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
#' @param genes Optional character vector of gene symbols to restrict the
#'   download (e.g., \code{c("TP53", "BRCA1", "EGFR")}). Only used when
#'   \code{source = "cBioPortal"}. If \code{NULL} (default), all available
#'   genes are downloaded, which can be very slow for large studies.
#' @param custom_data Optional. A named list of matrices for custom data.
#'   Names must match \code{data_types}. Only used when \code{source = "custom"}.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}.
#'
#' @return A \code{MultiAssayExperiment} object containing the requested omics
#'   layers and associated clinical data.
#'
#' @examples
#' \dontrun{
#' # Load breast cancer data from cBioPortal (all genes)
#' brca <- load_cancer_data("brca_tcga")
#'
#' # Load only selected genes (faster)
#' brca <- load_cancer_data(
#'   "brca_tcga",
#'   data_types = c("methylation", "expression"),
#'   genes = c("TP53", "BRCA1", "BRCA2", "EGFR")
#' )
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
                             genes = NULL,
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

    # Create cBioPortal API connection object
    cbio <- cBioPortalData::cBioPortal()

    # --- Discover available molecular profiles for this study ---
    if (verbose) message("  Discovering available molecular profiles...")

    all_profiles <- tryCatch(
      cBioPortalData::molecularProfiles(cbio, cancer_type),
      error = function(e) {
        stop("Could not retrieve profiles for study '", cancer_type,
             "'. Check that the study ID is correct. Error: ",
             conditionMessage(e))
      }
    )

    profile_ids_all <- all_profiles$molecularProfileId

    # --- Find matching profile IDs for requested data types ---
    profile_ids <- character(0)
    profile_names <- character(0)

    if ("methylation" %in% data_types) {
      hit <- profile_ids_all[grepl("methylation", profile_ids_all,
                                   ignore.case = TRUE)]
      hit <- hit[!grepl("median|zscore|zscores", hit, ignore.case = TRUE)]
      if (length(hit) > 0) {
        profile_ids <- c(profile_ids, hit[1])
        profile_names <- c(profile_names, "methylation")
      } else {
        warning("No methylation profile found for study '", cancer_type, "'.")
      }
    }

    if ("expression" %in% data_types) {
      hit <- profile_ids_all[grepl("rna_seq_v2_mrna$", profile_ids_all,
                                   ignore.case = TRUE)]
      if (length(hit) == 0) {
        hit <- profile_ids_all[grepl("_mrna$", profile_ids_all,
                                     ignore.case = TRUE) &
                                 !grepl("median|zscore|zscores",
                                        profile_ids_all, ignore.case = TRUE)]
      }
      if (length(hit) > 0) {
        profile_ids <- c(profile_ids, hit[1])
        profile_names <- c(profile_names, "expression")
      } else {
        warning("No expression profile found for study '", cancer_type, "'.")
      }
    }

    if ("mirna" %in% data_types) {
      hit <- profile_ids_all[grepl("mirna", profile_ids_all,
                                   ignore.case = TRUE)]
      hit <- hit[!grepl("median|zscore|zscores", hit, ignore.case = TRUE)]
      if (length(hit) > 0) {
        profile_ids <- c(profile_ids, hit[1])
        profile_names <- c(profile_names, "mirna")
      } else {
        warning("No miRNA profile found for study '", cancer_type, "'.")
      }
    }

    if (length(profile_ids) == 0) {
      stop("None of the requested data types (",
           paste(data_types, collapse = ", "),
           ") are available for study '", cancer_type, "'.")
    }

    if (verbose) {
      message("  Using profiles:")
      for (i in seq_along(profile_ids)) {
        message("    ", profile_names[i], " -> ", profile_ids[i])
      }
    }

    # --- Determine gene list ---
    if (is.null(genes)) {
      # Download all available genes (may be slow)
      if (verbose) message("  Retrieving full gene list (this may take a moment)...")
      gene_table <- tryCatch(
        cBioPortalData::geneTable(cbio),
        error = function(e) {
          stop("Could not retrieve gene table: ", conditionMessage(e))
        }
      )
      gene_symbols <- unique(gene_table$hugoGeneSymbol)
      if (verbose) message("  Found ", length(gene_symbols), " genes.")
    } else {
      gene_symbols <- genes
      if (verbose) message("  Using user-provided gene list (",
                           length(gene_symbols), " genes).")
    }

    # --- Load the data ---
    mae <- tryCatch(
      cBioPortalData::cBioPortalData(
        api = cbio,
        studyId = cancer_type,
        genes = gene_symbols,
        by = "hugoGeneSymbol",
        molecularProfileIds = profile_ids
      ),
      error = function(e) {
        stop("Failed to load data from cBioPortal: ", conditionMessage(e))
      }
    )

    # --- Rename assays to user-friendly names ---
    if (length(profile_names) ==
        length(names(MultiAssayExperiment::experiments(mae)))) {
      names(MultiAssayExperiment::experiments(mae)) <- profile_names
    }

    if (verbose) message("Successfully loaded ",
                         length(MultiAssayExperiment::experiments(mae)),
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
