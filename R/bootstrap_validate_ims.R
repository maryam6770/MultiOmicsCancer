#' Bootstrap validation of an Immune Methylation Score
#'
#' Estimates the optimism of an IMS built with \code{\link{build_ims}}
#' using Harrell's bootstrap method. The entire pipeline (M-value
#' transformation, variance filter, Elastic Net fitting) is repeated in
#' each bootstrap replicate, and the optimism-corrected C-index is
#' reported.
#'
#' @param mae A \code{MultiAssayExperiment} or \code{SummarizedExperiment}.
#' @param n_boot Number of bootstrap replicates (default 200).
#' @param verbose Logical; print progress.
#' @param ... Additional arguments passed to \code{\link{build_ims}}.
#'
#' @return A \code{BootstrapIMSResult} object with elements
#'   \code{c_apparent}, \code{optimism}, \code{c_corrected},
#'   \code{c_boot}, \code{n_boot}, and \code{parameters}.
#'
#' @examples
#' data(example_mae)
#' set.seed(42)
#' bv <- bootstrap_validate_ims(
#'     example_mae,
#'     immune_genes = c("CD274", "PDCD1", "BTLA", "TNF", "CD8A"),
#'     n_boot = 20,
#'     verbose = FALSE
#' )
#' bv
#' @export
bootstrap_validate_ims <- function(mae,
                                    n_boot = 200L,
                                    verbose = TRUE,
                                    ...) {
    args <- list(...)
    time_column <- if (!is.null(args$time_column)) args$time_column else "OS.time"
    event_column <- if (!is.null(args$event_column)) args$event_column else "OS"
    assay_name <- if (!is.null(args$assay_name)) args$assay_name else "methylation"
    immune_genes <- if (!is.null(args$immune_genes)) args$immune_genes else .default_immune_genes()
    gene_map <- if (!is.null(args$gene_map)) args$gene_map else NULL
    min_sd <- if (!is.null(args$min_sd)) args$min_sd else 0.1
    alpha <- if (!is.null(args$alpha)) args$alpha else 0.5
    nfolds <- if (!is.null(args$nfolds)) args$nfolds else 5L

    clin <- as.data.frame(MultiAssayExperiment::colData(mae))
    if (!time_column %in% colnames(clin)) {
        stop("Column '", time_column, "' not found in colData.", call. = FALSE)
    }
    if (!event_column %in% colnames(clin)) {
        stop("Column '", event_column, "' not found in colData.", call. = FALSE)
    }
    time <- as.numeric(clin[[time_column]])
    event <- as.numeric(clin[[event_column]])

    se <- .extract_assay(mae, assay_name)
    beta_all <- SummarizedExperiment::assay(se, assay_name)

    if (is.null(gene_map)) {
        rd <- SummarizedExperiment::rowData(se)
        if (!"gene" %in% colnames(rd)) {
            stop("No 'gene' column in rowData and no gene_map supplied.", call. = FALSE)
        }
        gene_map <- data.frame(cpg = rownames(se), gene = as.character(rd[["gene"]]),
                               stringsAsFactors = FALSE)
    } else if (is.character(gene_map) && !is.null(names(gene_map))) {
        gene_map <- data.frame(cpg = names(gene_map), gene = unname(gene_map),
                               stringsAsFactors = FALSE)
    }
    keep_gene <- gene_map$gene %in% immune_genes
    gene_map <- gene_map[keep_gene, , drop = FALSE]
    beta_all <- beta_all[gene_map$cpg, , drop = FALSE]

    ref <- do.call(build_ims, c(list(mae = mae, verbose = FALSE), args))
    c_apparent <- .ims_cindex(ref$score, time, event)

    if (verbose) message(sprintf("Apparent C-index: %.3f", c_apparent))

    n <- nrow(clin)
    c_boot <- rep(NA_real_, n_boot)
    optimism <- rep(NA_real_, n_boot)

    for (b in seq_len(n_boot)) {
        idx <- sample.int(n, n, replace = TRUE)
        res_b <- tryCatch(
            .ims_pipeline(beta_all[, idx, drop = FALSE], gene_map,
                           time[idx], event[idx],
                           min_sd = min_sd, alpha = alpha, nfolds = nfolds),
            error = function(e) NULL
        )
        if (is.null(res_b)) next
        c_boot[b] <- .ims_cindex(res_b$score, time[idx], event[idx])
        score_ref <- .ims_apply_ref(ref, beta_all[, idx, drop = FALSE])
        c_orig_b <- .ims_cindex(score_ref, time[idx], event[idx])
        optimism[b] <- c_boot[b] - c_orig_b
        if (verbose && b %% 25 == 0) message(sprintf("  bootstrap %d/%d", b, n_boot))
    }

    c_corrected <- c_apparent - mean(optimism, na.rm = TRUE)

    out <- list(
        c_apparent = c_apparent,
        optimism = mean(optimism, na.rm = TRUE),
        c_corrected = c_corrected,
        c_boot = c_boot,
        n_boot = n_boot,
        parameters = list(
            time_column = time_column, event_column = event_column,
            n_samples = n, n_events = sum(event, na.rm = TRUE)
        )
    )
    class(out) <- c("BootstrapIMSResult", "list")
    out
}

#' @export
print.BootstrapIMSResult <- function(x, ...) {
    cat("BootstrapIMSResult\n")
    cat("  Apparent C-index:   ", round(x$c_apparent, 3), "\n")
    cat("  Mean optimism:      ", round(x$optimism, 3), "\n")
    cat("  Corrected C-index:  ", round(x$c_corrected, 3), "\n")
    cat("  Bootstrap reps:     ", x$n_boot, "\n")
    invisible(x)
}

.ims_cindex <- function(score, time, event) {
    ok <- !is.na(score) & !is.na(time) & !is.na(event) & time > 0
    if (sum(ok) < 5 || sum(event[ok]) < 3) return(NA_real_)
    fit <- survival::coxph(survival::Surv(time[ok], event[ok]) ~ score[ok])
    as.numeric(summary(fit)$concordance[1])
}

.ims_pipeline <- function(beta, gene_map, time, event, min_sd, alpha, nfolds) {
    beta[beta <= 0] <- 0.001; beta[beta >= 1] <- 0.999
    if (anyNA(beta)) {
        for (j in seq_len(ncol(beta))) {
            if (anyNA(beta[, j])) beta[is.na(beta[, j]), j] <- stats::median(beta[, j], na.rm = TRUE)
        }
    }
    mval <- log2(beta / (1 - beta))
    sd_m <- apply(mval, 1, stats::sd)
    keep <- !is.na(sd_m) & sd_m > min_sd
    if (sum(keep) < 3) return(NULL)
    mval <- mval[keep, , drop = FALSE]
    gm <- gene_map[keep, , drop = FALSE]
    X <- t(mval)
    Y <- survival::Surv(time, event)
    ok <- !is.na(time) & !is.na(event) & time > 0
    X <- X[ok, , drop = FALSE]; Y <- Y[ok]
    foldid <- sample(rep(seq_len(nfolds), length.out = nrow(X)))
    fit <- tryCatch(.cv_glmnet_safe(X, Y, family = "cox", cox.ties = "breslow", alpha = alpha, foldid = foldid),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    coefs <- as.numeric(stats::coef(fit, s = "lambda.min"))
    names(coefs) <- colnames(X)
    if (sum(coefs != 0) == 0) return(NULL)
    score <- as.numeric(X %*% coefs)
    names(score) <- rownames(X)
    list(score = score, coefficients = coefs[coefs != 0], gene_map = gm)
}

.ims_apply_ref <- function(ref, beta) {
    cpgs <- names(ref$coefficients)
    cpgs <- intersect(cpgs, rownames(beta))
    if (length(cpgs) == 0) return(rep(NA_real_, ncol(beta)))
    b <- beta[cpgs, , drop = FALSE]
    b[b <= 0] <- 0.001; b[b >= 1] <- 0.999
    m <- log2(b / (1 - b))
    m[is.na(m)] <- 0
    as.numeric(ref$coefficients[cpgs] %*% m)
}
