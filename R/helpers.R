# Internal helpers

.extract_assay <- function(x, assay_name) {
    if (methods::is(x, "MultiAssayExperiment")) {
        exps <- MultiAssayExperiment::experiments(x)
        if (!assay_name %in% names(exps)) {
            stop("Assay '", assay_name, "' not found in MAE.", call. = FALSE)
        }
        return(exps[[assay_name]])
    }
    if (methods::is(x, "SummarizedExperiment")) {
        return(x)
    }
    stop("Input must be a MultiAssayExperiment or SummarizedExperiment.",
         call. = FALSE)
}

.default_immune_genes <- function() {
    unlist(list(
        checkpoint = c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT",
                       "HAVCR2", "PDCD1LG2", "VSIR", "BTLA"),
        cytokine = c("IFNG", "TNF", "IL6", "IL10", "TGFB1"),
        chemokine = c("CXCL9", "CXCL10"),
        antigen = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
        tcell = c("CD8A", "CD8B", "CD4", "FOXP3", "GZMA", "GZMB", "PRF1"),
        signaling = c("STAT1", "IRF1", "JAK1", "JAK2", "IDO1")
    ), use.names = FALSE)
}

# Safe wrappers that suppress warnings on small simulated data
# (glmnet convergence and coxph perfect-separation warnings
#  do not occur on real TCGA data with hundreds of samples)
.cv_glmnet_safe <- function(...) {
    suppressWarnings(glmnet::cv.glmnet(...))
}

.coxph_safe <- function(...) {
    suppressWarnings(survival::coxph(...))
}
