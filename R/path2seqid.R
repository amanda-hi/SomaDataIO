#' Convert Pathway Term to SeqIds
#'
#' Utility function to convert a GO term or other biological/functional pathway
#' to a vector of SeqIds associated with the features (most commonly genes) in
#' that pathway.
#'
#' @param path Character. A pathway identifier or vector of pathway identifiers.
#' Pathway identifiers must be from the GO or MSigDB data resources.
#' @param format Character. The database to which the value in `path` belongs.
#' Accepted values are "GO" or "MSigDB" (case insensitive).
#' @return A list containing a vector of SeqIds. If multiple identifiers are
#' provided to `path`, the returned list will contain a sub-list for each
#' identifier.
#' @author Amanda Hiser
#' @examples
#' # Single GO term
#' path2seqid("GO:0007596")
#'
#' # Vector of GO terms
#' path2seqid(c("GO:0007596", "GO:0098927", "GO:0002065"))
#'
#' # MSigDB gene set
#' path2seqid("M2622_GLI1_UP.V1_DN", format = "msigdb")
#'
#' # Pathways are case-insensitive and can be used as keywords for a search
#' path2seqid("myeloid", format = "msigdb") |> head(n = 2)
#' @export
path2seqid <- function(path, format = c("go", "msigdb")) {

  format <- match.arg(format)

  if ( format == "go" ) {
    # Downside is it loses information about which category it comes from
    all_go <- c(go2seqid$biological_process, go2seqid$molecular_function)

    # Strip prefixes & other characters
    pattern <- gsub("\\D", "", path)
    idx <- sapply(pattern, function(x) grep(x, names(all_go)))

    if ( length(idx) == 0L ) {
      stop(
        paste("Couldn't find SeqIds associated with the supplied GO term(s):", path,
              "\nIt may be misspelled, or may belong to the 'Cellular Component'",
              "category, which is currently not supported."),
        call. = FALSE) # Could use stopifnot length >= 1, but this is more informative
    } else {
      match <- all_go[unname(idx)]
    }
  }

  if ( format == "msigdb" ) {
    all_msig <- unlist(msig2seqid, recursive = FALSE)
    pattern <- toupper(path) # All MSigDB gene set names are in uppercase
    idx <- sapply(pattern, function(x) grep(x, names(all_msig)))
    idx <- unlist(idx)

    stopifnot("Couldn't find SeqIds associated with the supplied MSigDB gene set(s)" = length(idx) >= 1L)

    match <- all_msig[unname(idx)]
  }

  match
}
