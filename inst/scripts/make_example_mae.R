# make_example_mae.R
# Generates data/example_mae.rda used in the MultiOmicsCancer vignette
# and examples. Not run during package build; kept for reproducibility.

set.seed(123)

library(MultiAssayExperiment)
library(SummarizedExperiment)

n_tumor  <- 20
n_normal <- 10
n_total  <- n_tumor + n_normal

sample_ids <- c(paste0("TUMOR_", 1:n_tumor), paste0("NORMAL_", 1:n_normal))
patient_ids <- sample_ids
sample_type <- c(rep("Primary Tumor", n_tumor), rep("Solid Tissue Normal", n_normal))

genes <- c("TP53", "BRCA1", "EGFR", "CD274", "PDCD1", "BTLA", "TNF", "CD8A")
cpgs  <- paste0("cg", sprintf("%08d", 1:50))

expr_mat <- matrix(rnorm(length(genes) * n_total, mean = 5, sd = 2),
                   nrow = length(genes), ncol = n_total)
rownames(expr_mat) <- genes
colnames(expr_mat) <- sample_ids

meth_mat <- matrix(runif(length(cpgs) * n_total, 0.1, 0.9),
                   nrow = length(cpgs), ncol = n_total)
rownames(meth_mat) <- cpgs
colnames(meth_mat) <- sample_ids
meth_mat[1:10, 1:n_tumor] <- pmin(pmax(meth_mat[1:10, 1:n_tumor] + 0.2, 0.01), 0.99)

clin_df <- data.frame(
    sampleID    = sample_ids,
    patientID   = patient_ids,
    sample_type = sample_type,
    OS.time     = round(runif(n_total, 100, 2000)),
    OS          = sample(c(0, 1), n_total, replace = TRUE, prob = c(0.7, 0.3)),
    stringsAsFactors = FALSE
)
rownames(clin_df) <- sample_ids

se_expr <- SummarizedExperiment(
    assays  = list(expression = expr_mat),
    colData = DataFrame(sampleID = sample_ids, sample_type = sample_type)
)
se_meth <- SummarizedExperiment(
    assays  = list(methylation = meth_mat),
    colData = DataFrame(sampleID = sample_ids, sample_type = sample_type)
)

example_mae <- MultiAssayExperiment(
    experiments = list(expression = se_expr, methylation = se_meth),
    colData     = DataFrame(clin_df)
)

dir.create("data", showWarnings = FALSE)
save(example_mae, file = "data/example_mae.rda", compress = "xz")
