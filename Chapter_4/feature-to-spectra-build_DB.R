## -------------------------------------------------------------------------------------------
## Fuctions used for feature-to-spectra matching
## -------------------------------------------------------------------------------------------

# build_DB ----------------------------------------------------------------
## function for building a database table from RDA files
build_DB <- function(wd, # path to RDA files
                     spec.threshold, # BPI threshold for spectrum features
                     chem.threshold # BPI threshold for spectrum features
) {
  # set default intensity threshold spectrum-wise
  if (missing(spec.threshold)) {
    spec.threshold <- 0
    message("spec.threshold value not selected, using default: 0")
  }
  # set default intensity threshold chemical-wise
  if (missing(chem.threshold)) {
    chem.threshold <- 0
    message("chem.threshold value not selected, using default: 0")
  }

  # list DB compound rda files
  file.list <- list.files(path = wd, pattern = "*.rda", full.names = T)

  all.chem <- NULL

  # get progress bar
  pb <- txtProgressBar(min = 0, max = length(file.list), style = 3)


  #### for each DB compound -----
  for (i in 1:length(file.list)) {
    # determine the number of spectra existing in the chem.file
    load(file = file.list[i])
    spectra <- length(chem.file$analytical)
    all.spec <- NULL

    #### for each spectrum of the compound ----
    for (j in 1:spectra) {
      # normalise all features IT within spectrum to %BPI and remove features with IT < spec.threshold
      spec.norm.it <- ((chem.file$analytical[[j]]$spec[, "IT"]) / max(chem.file$analytical[[j]]$spec[, "IT"])) * 100
      spec <- cbind(chem.file$analytical[[j]]$spec, spec.norm.it)
      spec <- spec[(which(spec[, "spec.norm.it"] > spec.threshold)), ]

      # create vector of repeated RT vals for each IT passing feature in the spec
      spec.rt <- rep((chem.file$analytical[[j]]$rt), nrow(spec))
      # create vector of repeated spectrum # vals for each IT passing feature in the spec
      spec.no <- rep(j, nrow(spec))
      # bind the passing spectrum features with the RT and spec # vals
      spec.final <- cbind(spec, spec.rt, spec.no)
      # creates a vector of unique peaks in the spectrum
      spec.final$peak.no <- 1:nrow(spec.final)
      # appends spec.final matrix to the matrix for all spectra
      all.spec <- rbind(all.spec, spec.final)
    }

    #### normalise all feature IT within chemical (all spectra) to %BPI and remove spectrum features with IT < chem.threshold ----
    chem.norm.it <- ((all.spec[, "IT"]) / max(all.spec[, "IT"])) * 100
    chem <- cbind(all.spec, chem.norm.it)
    chem <- chem[(which(chem[, "chem.norm.it"] > chem.threshold)), ]

    #### find major rt group ----
    # find sum of norm.it for each rt group
    chem.norm.it.t <- unlist(lapply(unique(chem$spec.no), function(s) {
      sum(subset(chem, spec.no == s, select = chem.norm.it))
    }))
    major.spec <- unique(chem$spec.no)[which.max(chem.norm.it.t)]
    chem$major.spec <- rep(FALSE, nrow(chem))
    chem[which(chem$spec.no == major.spec), "major.spec"] <- rep(TRUE, length(which(chem$spec.no == major.spec)))

    #### add chemical DB  details ----
    chem.no <- rep(i, nrow(chem))
    chem.NPCID <- rep(chem.file$id$NPCID, nrow(chem))
    chem.name <- rep(chem.file$id$name, nrow(chem))

    chem.final <- cbind(chem, chem.no, chem.NPCID, chem.name)
    all.chem <- rbind(all.chem, chem.final)

    #### update progress bar ----
    setTxtProgressBar(pb, i)
  }

  row.names(all.chem) <- NULL
  return(all.chem)
  close(pb)
}
