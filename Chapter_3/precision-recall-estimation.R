## -------------------------------------------------------------------------------------------
## Fuctions used for precision/recall estimation
## -------------------------------------------------------------------------------------------

# prep_xcms ---------------------------------------------------------------
## prepare final consensus map obtained by the XCMS
## function takes xcmsSet object as an argument
prep_xcms <- function(object) {
  sets_ind <- object@groupidx
  sets <- lapply(seq_len(length(sets_ind)), function(j) {
    set_ind <- sets_ind[[j]]
    set_samples <- data.frame(
      sample = object@peaks[set_ind, "sample"],
      peakid = object@peaks[set_ind, "peakid"],
      row.names = NULL
    )
    set_peakids <- unique(object@peaks[set_ind, c("peakid")])
    return(list(set_peakids, set_samples))
  })
  return(sets)
}

# prep_massFlowR ----------------------------------------------------------
## prepare final consensus map obtained by massFlowR
## function takes massFlowTemplate class object
prep_massFlowR <- function(object, valid = FALSE, peakid_original = FALSE) {
  ## if peaks were removed during data simulation, use column 'peakid_original' instead of the default 'peakid'
  peakid_colname <- ifelse(peakid_original, "peakid_original", "peakid")
  sets_ind <- if (valid == TRUE) {
    object@valid$peakid
  } else {
    object@tmp$peakid
  }
  lapply(seq(length(sets_ind)), function(j) {
    set_ind <- sets_ind[[j]]
    set_samples <- lapply(1:length(object@data), function(s) {
      if (set_ind %in% object@data[[s]]$tmp_peakid) {
        data.frame(
          sample = s,
          peakid = object@data[[s]][[peakid_colname]][match(set_ind, object@data[[s]]$tmp_peakid)],
          row.names = NULL
        )
      }
    })
    set_samples <- do.call("rbind", set_samples)
    return(list(unique(set_samples$peakid), set_samples))
  })
}

# prep_ground -------------------------------------------------------------
## prepare ground truth - the optimal consensus map
prep_ground <- function(sample_map, maps) {
  ## ground truth is a single list with N entries (N - number of consensus features)
  ## each list has S entries (S - number of maps in the experiment)
  ## consensus features are extracted from the sample map (which can be the first map in the experiment, provided that it contains all consensus features)
  cfeatures <- sample_map$peakid
  lapply(cfeatures, function(i) {
    lapply(seq_len(length(maps)), function(m) {
      map <- maps[[m]]
      i %in% map$peakid
    })
  })
}

# eval --------------------------------------------------------------------
## evaluate the alignment results against the ground truth
eval <- function(consensus_map, ground_truth) {
  precision <- 0
  recall <- 0

  ## for every consensus feature in the ground truth
  N <- length(ground_truth)
  for (i in seq_len(N)) {
    # print(i)
    ## |gt_i|, the number of features corresponding to the consensus feature in the ground truth
    gt_i_len <- length(which(unlist(ground_truth[[i]])))

    ## find which feature sets obtained by the tool contain the consensus feature
    tool_ind <- which(sapply(consensus_map, function(j) {
      i %in% j[[1]]
    }))

    if (length(tool_ind) > 0) {
      tool_i <- consensus_map[tool_ind]
      ## |tool_i|, number of features in the same group as the consensus features
      tool_i_len <- sum(sapply(tool_i, function(s) {
        nrow(s[[2]])
      }))

      ## |gt_i union tool_i|, the number of features that correspond to the consensus feature in consensus map
      gt_tool_i_len <- sum(sapply(tool_i, function(ij) {
        length(which(ij[[2]][, 2] == i))
      }))
      ## the number of features groups with the consensus feature
      m_i_len <- length(tool_ind)

      precision_i <- gt_tool_i_len / tool_i_len
      recall_i <- gt_tool_i_len / m_i_len / gt_i_len

      precision <- precision + precision_i
      recall <- recall + recall_i
    }
  }
  return(list(
    precision = precision / N,
    recall = recall / N
  ))
}
