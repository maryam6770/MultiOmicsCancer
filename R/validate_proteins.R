#' Validate Genes at Protein Level
#'
#' @description
#' Validates genes using UniProt (protein existence).
#'
#' @param genes Character vector of gene symbols.
#' @param use_uniprot Logical. Query UniProt? Default TRUE.
#' @param timeout Numeric. HTTP timeout in seconds. Default 30.
#'
#' @return A data frame with columns: gene, protein_exists, protein_name.
#'
#' @examples
#' # Source: Mock data for offline demonstration.
#' genes <- c("BRCA1", "TP53", "CD274")
#' res <- validate_proteins(genes, use_uniprot = FALSE)
#' head(res)
#' @export
validate_proteins <- function(genes, use_uniprot = TRUE, timeout = 30) {
  if (!is.character(genes) || length(genes) == 0) {
    stop('genes must be a non-empty character vector.')
  }
  res <- data.frame(
    gene = genes,
    protein_exists = NA,
    protein_name = NA_character_,
    stringsAsFactors = FALSE
  )
  if (use_uniprot) {
    for (i in seq_along(genes)) {
      g <- genes[i]
      url <- paste0(
        'https://rest.uniprot.org/uniprotkb/search?',
        'query=gene_exact:', g, '%20AND%20organism_id:9606',
        '&format=json&size=1'
      )
      out <- tryCatch({
        resp <- httr2::request(url) |>
          httr2::req_timeout(timeout) |>
          httr2::req_perform()
        json <- httr2::resp_body_json(resp)
        if (length(json[['results']]) > 0) {
          entry <- json[['results']][[1]]
          pname <- tryCatch({
            v <- entry[['proteinDescription']][['recommendedName']][['fullName']][['value']]
            if (is.null(v)) NA_character_ else as.character(v)
          }, error = function(e) NA_character_)
          list(exists = TRUE, name = pname)
        } else {
          list(exists = FALSE, name = NA_character_)
        }
      }, error = function(e) list(exists = NA, name = NA_character_))
      res[['protein_exists']][i] <- out[['exists']]
      res[['protein_name']][i] <- out[['name']]
      Sys.sleep(0.2)
    }
  }
  res
}

