# To get things up and running for testing, I'm going to use the cleaned
# networks from SGRQ.
# NOTE: these networks were made using 11k

##############################
#--- Setup ------------------#
##############################

# Seems like the easiest way to do this (for now), although it
# introduces a dependency that could be avoided by downloading
# the database to an XML like SGRQ did
library(msigdb)

go2seqid <- list() # Gene Ontology (GO)
msig2seqid <- list() # Molecular Signatures Database (MSigDB)

checkForDups <- function(db) {
  res <- lapply(db, function(x) {
    sapply(x[[1]], duplicated)
  }) |> unlist()

  sum(res) == 0L
}


##############################
#--- GO Terms ---------------#
##############################

# Going to have to find a way to access this information without using SGRQ.
# Could possibly store it as a data object in SomaScan.db, then access it
# from there?
# Or re-generate it using the same method John did in SGRQ, and save this
# script for documentation and just load an object of all the pathways.

# Checking to see if these have been filtered for network size yet
# Verdict: they have not
# somaGetRichQuick::networks$GO_biological_process |>
#   group_by(go_bp_id) |>
#   count() |>
#   arrange(desc = n) |>
#   head()
#
# somaGetRichQuick::networks$GO_molecular_function |>
#   group_by(go_mf_id) |>
#   count() |>
#   arrange(desc = n) |>
#   head()
#
# somaGetRichQuick::networks$MSigDB_immunologic_signature |>
#   group_by(systematic_name) |>
#   count() |>
#   arrange(desc = n) |>
#   head()

go2seqid$biological_process <- somaGetRichQuick::networks$GO_biological_process
go2seqid$molecular_function <- somaGetRichQuick::networks$GO_molecular_function

# Create single "term" col and remove unneeded cols
go2seqid <- lapply(go2seqid, function(x) {
  id_col <- grep("\\_id", colnames(x))
  name_col <- grep("\\_name", colnames(x))
  x$term <- paste(x[,id_col], x[,name_col], sep = "_")
  x$term <- gsub("GO\\:", "", x$term) # Remove prefix to comply with R's var name rules
  x$term <- gsub("\\s+", "_", x$term)
  # !!!!! also need to replace other non-R characters, like / and - and ., otherwise the name gets wrapped in ticks anyway

  if ( any(grepl("\\s", x$term)) ) {
    stop("Spaces still present in GO term lists, exiting...")
  }

  x[,c("seqid", "term")]
})

# Iterate thru every GO term and subset original df to identify
# SeqIds associated with the term
go2seqid <- lapply(go2seqid, function(x) {
  v <- unique(x$term) |>
    setNames(paste0("GO_BP_", unique(x$term)))

  # Very slow, surely there is a better way to do this
  lapply(v, function(y) {
    x[x$term == y,]$seqid
  })
})

# Making sure none of the sets have duplicated values
stopifnot(checkForDups(go2seqid))


##############################
#--- MSigDb -----------------#
##############################

#----- From scratch ---------#

# Downloads the MSigDb, don't want to be running this often! Slow!
# test <- msigdb::getMsigdb(org = "hs", id = "EZID")
# listSubCollections(test)

# For licensing reasons, can't use collections containing any of the following:
# - KEGG
# - BioCarta
# - AAAS/STKE Cell Signaling Database
# This package already excludes KEGG, thankfully. So it could be a good
# starting point for aggregating gene sets, instead of using SGRQ to start

# These are found in the c2 collection
# cols <- listCollections(test)
# cols <- cols[cols != "c2"]
# test2 <- subsetCollection(test, collection = cols)
# listSubCollections(test2) # Making sure BioCarta is gone
# docs are here: https://bioconductor.org/packages/release/data/experiment/vignettes/msigdb/inst/doc/msigdb.html


#----- From SGRQ ---------#
# Will use this for now, just to get started. Eventually will be replaced
# by resource above (msigdb R package)
msigdb_idx <- grep("msigdb", names(somaGetRichQuick::networks), ignore.case = TRUE)
msig2seqid <- somaGetRichQuick::networks[msigdb_idx]

# Not using curated collection, contains KEGG & BioCarta
msig2seqid <- msig2seqid[-grep("curated", names(msig2seqid))]

# Create single "term" col and remove unneeded cols
msig2seqid <- lapply(msig2seqid, function(x) {
  x$term <- paste(x$systematic_name, x$standard_name, sep = "_")
  x$term <- gsub("\\s+", "_", x$term)

  if ( any(grepl("\\s", x$term)) ) {
    stop("Spaces still present in GO term lists, exiting...")
  }

  x[,c("seqid", "term")]
})

# Iterate thru every gene set and subset original df to identify SeqIds
# associated with the set
msig2seqid <- lapply(msig2seqid, function(x) {
  v <- unique(x$term) |>
    setNames(paste0("MSigDB_", unique(x$term)))

  # Very slow, surely there is a better way to do this
  lapply(v, function(y) {
    x[x$term == y,]$seqid
  })
})

# Making sure none of the sets have duplicated values
stopifnot(checkForDups(msig2seqid))


##############################
#--- Saving Objects ---------#
##############################

save(go2seqid, file = "data/go_seqIds.rda", compress = "xz")
save(msig2seqid, file = "data/msigdb_seqIds.rda", compress = "xz")
