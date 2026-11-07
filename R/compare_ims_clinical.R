#' Compare an IMS with a clinical model
#'
#' Compares the discrimination of a clinical model with and without the
#' IMS. Returns the ΔC-index, the likelihood-ratio test statistic, and
#' the two component C-indices.
#'
#' @param mae A \code{MultiAssayExperiment} or \code{SummarizedExperiment}.
#' @param ims_result An \code{IMSResult} from \code{\link{build_ims}}.
#' @param clinical_vars Character vector of clinical column names to use
#'   as covariates. If \code{NULL}, the IMS is compared against the
#'   intercept-only model.
#' @param time_column Column with survival time.
#' @param event_column Column with event indicator.
#' @param verbose Logical; print progress.
#'
#' @return A \code{IMSComparison} object with elements \code{c_clinical},
#'   \code{c_combined}, \code{delta_c}, \code{lrt_stat},
#'   \code{lrt_p}, \code{fit_clinical}, \code{fit_combined},
#'   \code{n_samples}, and \code{n_events}.
#'
#' @examples
#' data(example_mae)
#' set.seed(42)
#' ims <- build_ims(example_mae,
#'     immune_genes = c("CD274", "PDCD1", "BTLA", "TNF", "CD8A"),
#'     verbose = FALSE)
#' cmp <- compare_ims_clinical(example_mae, ims,
#'     clinical_vars = "sample_type",
#'     verbose = FALSE)
#' cmp
#' @export
compare_ims_clinical <- function(mae, ims_result,
                                  clinical_vars = NULL,
                                  time_column = "OS.time",
                                  event_column = "OS",
                                  verbose = TRUE) {
    if (!methods::is(ims_result, "IMSResult")) {
        stop("ims_result must be an IMSResult object.", call. = FALSE)
    }
    clin <- as.data.frame(MultiAssayExperiment::colData(mae))
    if (!time_column %in% colnames(clin)) {
        stop("Column '", time_column, "' not found in colData.", call. = FALSE)
    }
    if (!event_column %in% colnames(clin)) {
        stop("Column '", event_column, "' not found in colData.", call. = FALSE)
    }
    if (!is.null(clinical_vars)) {
        missing_vars <- setdiff(clinical_vars, colnames(clin))
        if (length(missing_vars) > 0) {
            stop("Clinical variables not found: ", paste(missing_vars, collapse = ", "),
                 call. = FALSE)
        }
    }

    score <- ims_result$score
    sids <- names(score)
    clin <- clin[sids, , drop = FALSE]
    time <- as.numeric(clin[[time_column]])
    event <- as.numeric(clin[[event_column]])
    ok <- !is.na(time) & !is.na(event) & time > 0 & !is.na(score)
    time <- time[ok]; event <- event[ok]; score <- score[ok]
    clin_ok <- clin[ok, , drop = FALSE]

    Y <- survival::Surv(time, event)

    if (is.null(clinical_vars) || length(clinical_vars) == 0) {
        fit_clin <- .coxph_safe(Y ~ 1)
    } else {
        f_clin <- stats::as.formula(paste("Y ~", paste(clinical_vars, collapse = " + ")))
        fit_clin <- .coxph_safe(f_clin, data = clin_ok)
    }

    f_comb <- stats::as.formula(paste("Y ~ score +",
        if (length(clinical_vars) > 0) paste(clinical_vars, collapse = " + ") else "1"))
    fit_comb <- .coxph_safe(f_comb, data = clin_ok)

    c_clin <- as.numeric(summary(fit_clin)$concordance[1])
    c_comb <- as.numeric(summary(fit_comb)$concordance[1])
    delta_c <- c_comb - c_clin

    lrt <- tryCatch(stats::anova(fit_clin, fit_comb, test = "LRT"),
                    error = function(e) NULL)
    if (!is.null(lrt) && nrow(lrt) >= 2) {
        lrt_stat <- as.numeric(lrt$Chisq[2])
        lrt_p    <- as.numeric(lrt$`Pr(>|Chi|)`[2])
    } else {
        lrt_stat <- NA_real_
        lrt_p    <- NA_real_
    }

    if (verbose) {
        message(sprintf("C (clinical): %.3f | C (clinical + IMS): %.3f | dC: %+.3f",
                        c_clin, c_comb, delta_c))
    }

    out <- list(
        c_clinical = c_clin,
        c_combined = c_comb,
        delta_c    = delta_c,
        lrt_stat   = lrt_stat,
        lrt_p      = lrt_p,
        fit_clinical = fit_clin,
        fit_combined = fit_comb,
        n_samples = length(time),
        n_events  = sum(event)
    )
    class(out) <- c("IMSComparison", "list")
    out
}

#' @export
print.IMSComparison <- function(x, ...) {
    cat("IMSComparison\n")
    cat("  C-index (clinical):          ", round(x$c_clinical, 3), "\n")
    cat("  C-index (clinical + IMS):    ", round(x$c_combined, 3), "\n")
    cat("  Delta C-index:               ", sprintf("%+.3f", x$delta_c), "\n")
    cat("  LRT statistic:               ", round(x$lrt_stat, 3), "\n")
    cat("  LRT p-value:                 ", format.pval(x$lrt_p, digits = 3), "\n")
    cat("  Samples / events:            ", x$n_samples, "/", x$n_events, "\n")
    invisible(x)
}
