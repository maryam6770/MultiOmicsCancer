#' Load Multi-Omics Cancer Data
#'
#' Loads multi-omics data for any cancer type from cBioPortal, TCGA, or
#' user-provided datasets. The function returns a \code{MultiAssayExperiment}
#' object that can be used by all downstream analysis functions in the package.
#'
#' @param cancer_type Character. Cancer study identifier. For cBioPortal, use
#'   study IDs such as \code{"brca_tcga"}, \code{"luad_tcga"}, \code{"coadread_tcga"}.
#'   For a full list, see \code{cBioPortalData::getStudies()}.
#' @param data_types Character vector. Omics layers to load. Options include
#'   \code{"methylation"}, \code{"expression"}, and \code{"mirna"}. Default loads
#'   all three.
#' @param source Character. Data source. One of \code{"cBioPortal"} (default),
#'   \code{"custom"}, or \code{"TCGA"}.
#' @param custom_data Optional. A named list of matrices for custom data.
#'   Names must match \code{data_types}. Only used when \code{source = "custom"}.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}.
#'
#' @return A \code{MultiAssayExperiment} object containing the requested omics
#'   layers and associated clinical data.
#'
#' @examples
#' \dontrun{
#' # Load breast cancer data from cBioPortal
#' brca <- load_cancer_data("brca_tcga")
#'
#' # Load only methylation and expression
#' luad <- load_cancer_data("luad_tcga",
#'                          data_types = c("methylation", "expression"))
#'
#' # Load custom data
#' custom <- load_cancer_data("my_study",
#'                            source = "custom",
#'                            custom_data = list(methylation = meth_matrix,
#'                                               expression = expr_matrix))
#' }
#'
#' @export
load_cancer_data <- function(cancer_type,
                             data_types = c("methylation", "expression", "mirna"),
                             source = "cBioPortal",
                             custom_data = NULL,
                             verbose = TRUE) {

  # --- Validate inputs ---
  if (!is.character(cancer_type) || length(cancer_type) != 1) {
    stop("'cancer_type' must be a single character string.")
  }

  valid_types <- c("methylation", "expression", "mirna")
  if (!all(data_types %in% valid_types)) {
    stop("'data_types' must be a subset of: ",
         paste(valid_types, collapse = ", "))
  }

  if (!source %in% c("cBioPortal", "custom", "TCGA")) {
    stop("'source' must be one of: 'cBioPortal', 'custom', 'TCGA'.")
  }

  # --- Case 1: User-provided custom data ---
  if (source == "custom") {
    if (is.null(custom_data)) {
      stop("When 'source = \"custom\"', 'custom_data' must be provided.")
    }
    if (!all(data_types %in% names(custom_data))) {
      stop("'custom_data' must contain all requested data types: ",
           paste(data_types, collapse = ", "))
    }

    if (verbose) message("Building MultiAssayExperiment from custom data...")

    experiments <- lapply(data_types, function(dt) {
      SummarizedExperiment::SummarizedExperiment(assays = list(counts = custom_data[[dt]]))
    })
    names(experiments) <- data_types

    mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = experiments)
    return(mae)
  }

  # --- Case 2: cBioPortal ---
  if (source == "cBioPortal") {
    if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
      stop("Package 'cBioPortalData' is required. Install with: ",
           "BiocManager::install('cBioPortalData')")
    }

    if (verbose) message("Loading data from cBioPortal for: ", cancer_type)

    # Map user-facing names to cBioPortal profile IDs
    profile_map <- c(
      methylation = "methylation",
      expression  = "mrna_seq_v2_rsem",
      mirna       = "mirna"
    )

    profile_ids <- unname(profile_map[data_types])

    mae <- tryCatch(
      cBioPortalData::cBioPortalData(
        studyId = cancer_type,
        molecularProfileIds = profile_ids
      ),
      error = function(e) {
        stop("Failed to load data from cBioPortal: ", conditionMessage(e))
      }
    )

    if (verbose) message("Successfully loaded ", length(experiments(mae)),
                         " omics layer(s).")
    return(mae)
  }

  # --- Case 3: TCGA via curatedTCGAData (optional dependency) ---
  if (source == "TCGA") {
    if (!requireNamespace("curatedTCGAData", quietly = TRUE)) {
      stop("Package 'curatedTCGAData' is required for TCGA data. ",
           "Install with: BiocManager::install('curatedTCGAData')")
    }

    if (verbose) message("Loading TCGA data for: ", cancer_type)

    mae <- curatedTCGAData::curatedTCGAData(
      diseaseCode = cancer_type,
      assays = data_types,
      dry.run = FALSE
    )
    return(mae)
  }
}
