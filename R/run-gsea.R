#' Perform Gene Set Enrichment Analysis (GSEA) with SomaScan
#'
#' @description
#'   Implements the algorithm from the `fgsea` package to perform gene set
#'   enrichment (GSEA), tailored specifically to SomaScan data.
#'
#' @section Min/Max Set Sizes:
#'   GSEA normalization is not accurate for very small feature sets - for
#'   example, for sets with >10 features, as few as 2-3 features can produce
#'   significant results. To avoid this, any feature sets smaller than
#'   `min_feats` will be discarded before performing GSEA. Please note that this threshold
#'   is applied _after_ filtering to remove features (i.e. SomaScan identifiers)
#'   from the feature set (gene list or other functional group) that are not
#'   present in the `ranks` data set. The default value (N = 15) is based on
#'   the minimum gene set threshold used in the Broad Institute's
#'   [GSEA software](https://www.gsea-msigdb.org/gsea/index.jsp).
#'
#' @param ranks A two column `data.frame` of statistical test results to be
#'   used as ranks. The data frame must include the following columns:
#'   \describe{
#'     \item{SeqId}{Character. SomaScan `SeqIds` or R-compatible
#'   `AptNames`. For more information about these identifier formats, please
#'   see `?SeqId`.}
#'     \item{Rank}{Numeric. Ranking values associated with each identifier.}
#'   }
#'   The suggested value for the `Rank` column is the absolute value of a test
#'   statistic, typically log2 fold-change values. However, _any_ prior
#'   statistical test result can be used, as long as the data frame contains
#'   SomaScan identifiers and a numeric metric that can be used for ranking.
#'
#'   Alternatively, a named numeric vector can also be used. The vector must
#'   contain ranking values associated with each vector. Names must be SomaScan
#'   identifiers.
#' @param pathIds List containing character vectors of SomaScan identifiers
#'   associated with pathways, gene sets, or functional groups of interest. The
#'   name of each list component should be the name of the pathway/gene set. The
#'   elements in each character vector should be SomaScan `SeqIds` or `AptNames`.
#'   See the [go2seqid] object for reference.
#' @param adat Original `soma_adat` object used to obtain statistical test
#'   results in `ranks`.
#' @param min_feats Numeric. Minimum number of features within a feature set
#'   (i.e. vectors of `SeqIds` in the list provided for `pathIds`) for that set
#'   to be retained. Default is 15. Sets smaller than this value will be discarded.
#'   If `NULL`, no minimum is used. See `Details` for more information.
#' @param max_feats Numeric. Maximum number of features within a category
#'   for that category to be retained. Default is `NULL`, meaning that no
#'   upper bound will be used. See `Details` for more information.
#' @param nperm Numeric. Number of permutations to perform. Default is 1000.
#' @return A `data.frame` object containing the results of gene set enrichment,
#'   including the following columns:
#' \item{pathway}{Name of the pathway, as in `names(pathIds)`}
#' \item{pval}{Enrichment p-value}
#' \item{padj}{FDR-adjusted p-value}
#' \item{log2err}{The expected error for the standard deviation of the p-value logarithm}
#' \item{ES}{Enrichment score, calculated via the [Broad GSEA implementation](https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Enrichment_Score_(ES))}
#' \item{NES}{Enrichment score normalized to mean enrichment of random samples of the same size}
#' \item{Size}{Size of the pathway _after_ removing `SeqIds` not present in `ranks` (this step is performed before GSEA)}
#' \item{leadingEdge}{Vector with indexes of leading edge genes that drive the enrichment, see [Running a Leading Edge Analysis](http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading)}
#' @author Amanda Hiser
#' @references A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert,
#'   M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, &
#'   J.P. Mesirov, Gene set enrichment analysis: A knowledge-based approach for
#'   interpreting genome-wide expression profiles,
#'   Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550,
#'   https://doi.org/10.1073/pnas.0506580102 (2005).
#' @examples
#' ranks <- t_tests[,c("SeqId", "rank")]
#' results <- run_gsea(ranks = ranks, pathIds = go2seqid$biological_process,
#'                     example_data)
#' head(results)
#' @importFrom fgsea fgseaMultilevel
#' @export
run_gsea <- function(ranks, pathIds, adat, min_feats = 15L, max_feats = NULL,
                     nperm = 1000) {

  assay_ver <- getSomaScanVersion(adat)
  checkSomaScanVersion(assay_ver)

  if ( is.data.frame(ranks) ) {
    stopifnot("`ranks` must be a 2-column data frame." = ncol(ranks) == 2L)

    # Identify column w/ SeqIds or AptNames
    seq_col <- which(
      apply(ranks, 2, function(x) all(is.apt(x)) | all(is.SeqId(x)))
    )
    ranks[[seq_col]] <- getSeqId(ranks[[seq_col]]) # Must be in SeqId format

    # Create ranks vector, if data frame is provided
    rank_vec <- ranks[[which(names(ranks) != names(seq_col))]]
    names(rank_vec) <- ranks[[seq_col]]

  } else if ( is.vector(ranks) ) {
    stopifnot("Names of `ranks` must be SomaLogic SeqIds or AptNames." = !is.apt(ranks) | !is.SeqId(ranks))

    rank_vec <- ranks
    names(rank_vec) <- getSeqId(names(ranks))
  }

  results <- fgsea::fgseaMultilevel(pathways = pathIds,
                                    stats = rank_vec,
                                    minSize = min_feats,
                                    maxSize = max_feats,
                                    nPermSimple = nperm,
                                    scoreType = "pos") # Reqd when all ranks >0 (abs value)

  results
}
