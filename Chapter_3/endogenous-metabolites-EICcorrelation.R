## -------------------------------------------------------------------------------------------
## Fuctions used for EIC correlation analysis
## -------------------------------------------------------------------------------------------

# Functions ---------------------------------------------------------------
prepROI <- function(data_dir, metabolites_ids = NULL, samples_fnames) {
  rtMin <- read.csv(list.files(data_dir, "rtMin", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  rtMax <- read.csv(list.files(data_dir, "rtMax", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  rt <- read.csv(list.files(data_dir, "rt.csv", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  mzMin <- read.csv(list.files(data_dir, "mzMin", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  mzMax <- read.csv(list.files(data_dir, "mzMax", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  mz <- read.csv(list.files(data_dir, "mz.csv", full.names = TRUE), header = TRUE, stringsAsFactors = FALSE, row.names = 1)

  tables <- list(rtMin, rtMax, rt, mzMin, mzMax, mz)

  ## extract specified metabolites using cpdID in the columns
  if (!is.null(metabolites_ids)) {
    tables <- lapply(tables, FUN = getROImetabolites, metabolites_ids = metabolites_ids)
  }

  ## extract specified samples because ppR output sample order can be different from originally supplied
  roi_all <- lapply(samples_fnames, FUN = getROI, tables = tables)
  return(roi_all)
}

getROImetabolites <- function(table, metabolites_ids) {
  as.data.frame(table[, match(metabolites_ids, colnames(table))])
}

getROI <- function(file_name, tables) {
  roi <- lapply(tables, function(table) {
    ind <- grep(file_name, rownames(table))
    if (length(ind) > 0) {
      as.data.frame(t(table[ind, ]))
    }
  })
  roi <- do.call("cbind", roi)
  if (!is.null(roi)) {
    roi <- setNames(roi,
      nm = c("rtMin", "rtMax", "rt", "mzMin", "mzMax", "mz")
    )
    rownames(roi) <- NULL
    return(roi)
  }
}
