load("data/go_msigdb_objects.rda")

checkForDups <- function(db) {
  res <- lapply(db, function(x) {
    sapply(x[[1]], duplicated)
  }) |> unlist()

  sum(res) == 0L
}

##############################
#--- GO Terms ---------------#
##############################

# Create single "term" col and remove unneeded cols
go2seqid <- lapply(go2seqid, function(x) {
  id_col <- grep("\\_id", colnames(x))
  name_col <- grep("\\_name", colnames(x))
  x$term <- paste(x[,id_col], x[,name_col], sep = "_")
  x$term <- gsub("GO\\:", "", x$term) # Remove prefix to comply with R's var name rules
  x$term <- gsub("\\s+", "_", x$term)
  # also need to replace other non-R characters, like / and - and ., otherwise
  # the name gets wrapped in ticks anyway

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
#--- Combined data frame ----#
##############################



##############################
#--- Saving Objects ---------#
##############################

save(go2seqid, file = "data/go_seqIds.rda", compress = "xz")
save(msig2seqid, file = "data/msigdb_seqIds.rda", compress = "xz")
