## -------------------------------------------------------------------------------------------
## Fuctions used for synthetic data generation
## -------------------------------------------------------------------------------------------

# simulateDATA ------------------------------------------------------------
## main function that simulates the desired number of peak-groups tables
simulateDATA <-
  function(fname,
           out_dir, # output directory
           meta_dir = NULL,
           mz_err = 0.001, # desired mz error
           rt_err = 2, # desired rt error
           drift_prop = 0, # desired proportion of peak-groups that experience RT drift
           drift = NULL, # path to .csv with nonlinear drift patters for reference compounds
           miss_probs = NULL, # path to .csv with probabilities of missingness for each peak
           miss = 0, # desired proportion of missing peaks
           batch_n = 0, # desired number of batches which resemble between-batch variation
           files_batch_n = 0 # desired number of simulated tables
  ) {
    if (!dir.exists(out_dir)) {
      ans <- 0
      while (ans < 1) {
        ans <- readline(
          paste0(
            "out_dir doesn't exist, create one? Enter Y/N \n"
          )
        )
        ## catch if input is N/n
        ans <- ifelse((grepl("N", ans) | grepl("n", ans)),
          2, 1
        )
        if (ans == 2) {
          stop("function was stopped.")
        } else {
          dir.create(out_dir, recursive = TRUE)
        }
      }
    }
    ## load the reference peak table that was acquired for a real QC sample
    ## the same peak table used for all simulations
    pks <- read.csv(fname, header = TRUE, stringsAsFactors = FALSE)

    ## calculate the number of peak-groups to be subjected to rt drift
    drift_n <- floor(length(unique(pks$peakgr)) * (drift_prop / 100))
    ## select peakgroups at random that will be subjected to rt drift
    ## for now, none of them are assigned for RT drift
    drift_peakgrs <-
      data.frame(peakgr = sort(sample(
        unique(pks$peakgr),
        size = drift_n, replace = F
      ), decreasing = F))

    ## should RT drift be applied to synthetic data?
    if (!is.null(drift)) {
      ## to model nonlinear drift using splines obtained on reference compounds
      drift <- read.csv(drift, header = TRUE, stringsAsFactors = FALSE)
      ## assign peakgroups to one of the drift patterns
      rcompounds <- unique(drift$rc)
      drift_peakgrs$drift <- sample(x = 1:length(rcompounds), size = drift_n, replace = TRUE)
      ## rt drift pattern is obtained by adding rt difference with the first rt in the splines
      drift_pattern <- lapply(rcompounds, function(rc) {
        drift_rc <- drift[which(drift$rc == rc), ]
        rtdiffs <- drift_rc$y - drift_rc$y[1]
        ## scale so that max to-first-sample difference would be equal to 4sec
        rtdiffs * 4 / max(abs(rtdiffs))
      })
    } else {
      ## if drift is not desired, apply natural log shape drift
      ## divide peakgroups into those that will elute earlier (- rt_err, 1), and later (+ rt_err, 2) in the run
      drift_peakgrs$drift <- sample(1:2, size = drift_n, replace = T)
      ## create rt drift pattern, which is added or substracted from the original rt of a peak-group
      drift_pattern <- lapply(c(-1, 1), function(d) {
        d * log(1:files_batch_n)
      })
    }

    ## should peaks be removed from the synthetic datasets?
    if (!is.null(miss_probs)) {
      ## load the peak missigness probabilities modelled earlier using modelMISS function
      probs <- read.csv(miss_probs, header = TRUE, stringsAsFactors = FALSE)
      ## how many peaks in the synthetic table should be missing, depends on user-selected proportion
      miss_n <- floor(nrow(pks) * miss)
    } else {
      probs <- NULL
      miss_n <- 0
    }
    ## should analytical batch structure be simulated?
    batches <- 1:ifelse(batch_n > 0, batch_n, 1)
    files_no <- seq(files_batch_n * max(batches))

    ## 1) Apply systematic rt drift to random peak-groups & random mz&rt variance to all peaks
    missed_peaks <- makeNOISE(
      pks = pks,
      out_dir = out_dir,
      mz_err = mz_err,
      rt_err = rt_err,
      drift_peakgrs = drift_peakgrs,
      drift_pattern = drift_pattern,
      probs = probs,
      miss_n = miss_n,
      batches = batches,
      files_batch_n = files_batch_n
    )

    # 2) Save metadata file for alignment
    filen <-
      sapply(1:files_batch_n, function(n) {
        paste0(paste0(rep(0, 3 - nchar(n)), collapse = ""), n)
      })
    batchn <- unlist(lapply(batches, function(n) {
      rep(n, files_batch_n)
    }))
    filenames <- paste0(
      "simulateDATA_batch", batchn, "_",
      filen
    )
    meta <- data.frame(
      filename = filenames,
      run_order = files_no,
      batch = batchn,
      raw_filepath = paste0(meta_dir, filenames, ".csv"),
      proc_filepath = paste0(meta_dir, filenames, ".csv"),
      stringsAsFactors = FALSE
    )
    write.csv(
      meta,
      file = paste0(out_dir, "metadata.csv"),
      quote = FALSE,
      row.names = FALSE
    )
    ## metadata local
    meta <- data.frame(
      filename = filenames,
      run_order = files_no,
      batch = batchn,
      raw_filepath = paste0(out_dir, filenames, ".csv"),
      proc_filepath = paste0(out_dir, filenames, ".csv"),
      stringsAsFactors = FALSE
    )
    write.csv(
      meta,
      file = paste0(out_dir, "metadata_local.csv"),
      quote = FALSE,
      row.names = FALSE
    )

    ## (3) Write which peak-groups were noisy
    if (nrow(missed_peaks) > 1) {
      drift_peakgrs <- merge(drift_peakgrs, missed_peaks)
    }
    write.csv(drift_peakgrs,
      file = paste0(out_dir, "noisy_peakgrs.csv"),
      quote = FALSE,
      row.names = FALSE
    )
    message("Tables were succesfully simulated")
  }

# makeNOISE ---------------------------------------------------------------
## simulates desired number of peak tables using provided RT drift and missingness patterns
makeNOISE <- function(pks,
                      out_dir,
                      mz_err,
                      rt_err,
                      drift_peakgrs,
                      drift_pattern,
                      probs,
                      miss_n,
                      batches,
                      files_batch_n) {
  nsamples <- length(batches) * files_batch_n

  ## generate intensities for n files, not for peaks marked for removal
  ints <- simulateINTENSITIES(pks = pks, nsamples = nsamples)

  ## create batch drift pattern
  # 1st batch is unchanged, other batches have either early or late elution in comparison to the 1st batch
  batch_gap <- c(0, sample(c(-10, 10), size = length(batches) - 1, replace = T))

  missed_peaks <- data.frame(
    peakgr = NA,
    peakid = NA,
    sample = NA,
    batch = NA
  )

  for (batch in batches) {
    message("Generating files for batch ", batch, "... ")

    ## a) Create the "template" for all files in this batch
    pks_batch <- pks
    ## introduce batch variance (if any)
    pks_batch$rt <- pks_batch$rt + batch_gap[batch]
    pks_batch$rtmin <- pks_batch$rtmin + batch_gap[batch]
    pks_batch$rtmax <- pks_batch$rtmax + batch_gap[batch]

    ## b) generate files_batch_n tables using the "template" peak table in this batch
    for (n in (1:files_batch_n)) {
      ## replace original intensity values with the simulated
      n_ori <- n + ifelse(batch == 1, 0, (batch - 1) * nsamples)
      pks_batch$into <- ints[, n_ori]

      ## introduce systematic rt variance to peak-groups (on top of batch effect)
      pks_n <- pks_batch
      for (peakgr in drift_peakgrs$peakgr) {
        peakgr_n <- which(drift_peakgrs$peakgr == peakgr)
        peakgr_rt <- pks_n[pks_n$peakgr == peakgr, c("rt", "rtmin", "rtmax")]
        ## what is the rt difference between the original rt and rt in this sample?
        peakgr_drift <- drift_pattern[[drift_peakgrs$drift[peakgr_n]]][n]
        pks_n[pks_n$peakgr == peakgr, "rt"] <- peakgr_rt$rt + peakgr_drift
        pks_n[pks_n$peakgr == peakgr, "rtmin"] <- peakgr_rt$rtmin + peakgr_drift
        pks_n[pks_n$peakgr == peakgr, "rtmax"] <- peakgr_rt$rtmax + peakgr_drift
      }
      ## introduce random mz&rt variance to all peaks
      pks_n <- makeRANDOMmz(pks = pks_n, mz_err = mz_err)
      pks_n <- makeRANDOMrt(pks = pks_n, rt_err = rt_err)

      ## remove peaks
      ## Probabilistic LOD, where the probability of a missing value increases at lower values
      if (!is.null(probs)) {
        P <- probs$P
        ## randomly select peaks to be removed depending on their missigness probability P
        missed <- sample(probs$peakid, size = miss_n, replace = FALSE, prob = P)
        missed <- probs[missed, ]

        ## remove peak-groups for which most intense peak was removed
        ## list either all peaks in a peakgroup (if most intense peak was randomly selected to be removed)
        ## or just the peaks that were randomly selected to be removed
        missed <- lapply(unique(missed$peakgr), function(pkg) {
          missed_pkg <- missed[missed$peakgr == pkg, ]
          probs_pks <- probs[probs$peakgr == pkg, ]
          peakids <- if (1 %in% missed_pkg$remove_order |
            nrow(missed_pkg) + 1 == nrow(probs_pks)) {
            list(probs$peakid[probs$peakgr == pkg], TRUE)
          } else {
            list(missed_pkg$peakid, FALSE)
          }
          data.frame(
            peakgr = pkg,
            peakid = peakids[[1]]
          )
        })
        ## save missing peaks id
        missed <- do.call("rbind", missed)
        missed$sample <- n
        missed$batch <- batch
        missed_peaks <- rbind(
          missed_peaks,
          missed
        )
        ## remove peaks
        pks_n <- pks_n[-missed$peakid, ]

        ## update peakids
        pks_n$peakid_original <- pks_n$peakid
        pks_n$peakid <- 1:nrow(pks_n)
      }

      ## write simulated table
      fname <- paste0(
        out_dir,
        "simulateDATA_batch",
        batch,
        "_",
        paste0(rep(0, 3 - nchar(n)), collapse = ""),
        n,
        ".csv"
      )
      write.csv(pks_n, file = fname, row.names = FALSE)
    }
  }
  return(missed_peaks)
}


# simulateINTENSITIES -----------------------------------------------------
## function to simulate intensities with pre-defined correlation structure
simulateINTENSITIES <- function(pks, pkgs = NULL, nsamples) {
  if (is.null(pkgs)) {
    ## for all peak groups in the sample table
    pkgs <- unique(pks$peakgr)
  }
  message("Simulating intensities for: ", length(pkgs), " peak-groups")

  ## add offset
  offset <- min(pks$into) + 1
  pks$into <- pks$into + offset

  ## generate a correlation matrix for every peak-group
  cor_list <- vector(mode = "list")
  for (pkg in pkgs) {
    cor_mat <- generateCORmat(pks = pks, pkg = pkg)
    cor_list[[length(cor_list) + 1]] <- cor_mat
  }
  ## create a correlation block matrix
  block_mat <- lapply(cor_list, "[[", "mat")
  block_mat <- as.matrix(Matrix::bdiag(block_mat))
  block_means <- unlist(lapply(cor_list, "[[", "means"))

  message("Finding nearest positive definite matrix ...")
  block_mat <- as.matrix(Matrix::nearPD(block_mat)$mat) # this take a while for all peaks

  message("Generating random variables ...")
  ## generate random variables
  dat <- mvtnorm::rmvnorm(n = nsamples, mean = block_means, sigma = block_mat, method = "chol")

  ## back transform
  dat <- exp(dat)

  dat <- as.data.frame(t(dat))
  colnames(dat) <- 1:nsamples
  dat$peakid <- unlist(lapply(cor_list, "[[", "peakids"))
  dat$peakgr <- rep(pkgs, unlist(lapply(cor_list, function(x) length(x[["peakids"]]))))

  ## arrange final table by the order in the original table
  dat <- dat[match(pks$peakid, dat$peakid), ]
  return(dat)
}

# generateCORmat --------------------------------------------------------------------------------------------------
## function to generate a correlation matrix for all peaks of a single peak group
generateCORmat <- function(pks, pkg) {
  peaks <- pks$peakid[which(pks$peakgr == pkg)]
  n <- length(peaks)
  ## diagonal correlations are very high
  cor_mat <- diag(sample(x = seq(0.9, 1 - sqrt(.Machine$double.eps), by = 0.00001), size = n), n, n)
  ## other correlations are lower
  cor_mat[upper.tri(cor_mat)] <- sample(x = seq(0.9, 1 - sqrt(.Machine$double.eps), by = 0.00001), size = (n * n - n) / 2)

  cor_mat <- Matrix::forceSymmetric(cor_mat)
  cor_mat <- as.matrix(Matrix::nearPD(cor_mat)$mat)
  means <- pks$into[peaks]

  ## convert correlations to covariance
  ## st devs are for arithmetic means scale
  stdevs <- means * 0.1
  stdevst <- stdevs %*% t(stdevs)
  cov_mat <- stdevst * cor_mat

  ## convert covariance matrix based on arithmetic returns to a covariance matrix from log-return
  cov_mat_log <- logreturn(means, cov_mat)
  means <- cov_mat_log$mean
  mat <- cov_mat_log$vcov
  return(list(peakids = peaks, means = means, mat = mat))
}

# logreturn ---------------------------------------------------------------
## code found @http://r.789695.n4.nabble.com/hep-on-arithmetic-covariance-conversion-to-log-covariance-td4646068.html
## function to convert arithmetic covariance to log-covariance
logreturn <- function(am, asigma) {
  M <- 1 / (1 + am)
  S <- log(diag(M) %*% asigma %*% diag(M) + 1)
  mu <- log(1 + am) - diag(S) / 2
  list(mean = mu, vcov = S)
}

# makeRANDOMmz ------------------------------------------------------------
## add random mz variance to all peaks, use user-selected mz variation range
makeRANDOMmz <- function(pks, mz_err) {
  pks_mz <- pks
  mz_noisy <- sapply(pks_mz$mz, function(x) {
    x + runif(min = -mz_err, max = mz_err, n = 1)
  })
  mz_dif <- pks_mz$mz - mz_noisy
  pks_mz$mz <- mz_noisy
  pks_mz$mzmin <- pks_mz$mzmin - mz_dif
  pks_mz$mzmax <- pks_mz$mzmax - mz_dif
  return(pks_mz)
}

# makeRANDOMrt ------------------------------------------------------------
## add random rt variance to peaks, use user-specified rt error
makeRANDOMrt <- function(pks, rt_err) {
  pks_rt <- pks
  rt_noisy <- sapply(pks_rt$rt, function(x) {
    x + runif(min = -rt_err, max = rt_err, n = 1)
  })
  rt_dif <- pks_rt$rt - rt_noisy
  pks_rt$rt <- rt_noisy
  pks_rt$rtmin <- pks_rt$rtmin - rt_dif
  pks_rt$rtmax <- pks_rt$rtmax - rt_dif
  return(pks_rt)
}

# modelMISS ---------------------------------------------------------------
## function to model missigness in the reference peak table
## the probability that a given peak will be missed in a simulated data inversely dependent on intensity
modelMISS <- function(fname, out_dir, miss) {
  ## load reference peak table
  pks <- read.csv(fname, header = TRUE, stringsAsFactors = FALSE)
  ints <- pks$into

  ## scale peak intensities from 0 to 1 gives the sharpest S shape of probability vs intensity
  ints <- scalePEAKS(ints, to = 1)

  ## given desired proportion of missing values, obtain optimal coefficient b0
  ## set coeficient beta1 to -10 to restrict the space for search
  b1 <- -10
  n <- length(ints)
  miss <- floor(miss * n)
  b0 <- optimiseBETA(n = n, b1 = b1, miss = miss, ints = ints)

  ## calculate missigness probabilities for every peak using its intensity values and obtained b0 and b1 values
  p_vals <- sapply(1:n, FUN = getPROB, ints = ints, b0 = b0, b1 = b1)

  ## obtain missigness probabilities, which are peak-group specific
  peaks <- lapply(unique(pks$peakgr), function(peakgr) {
    dt <- data.frame(
      peakgr = peakgr,
      peakid = pks$peakid[pks$peakgr == peakgr],
      ## peaks in a peak-group are removed according to their intensity (least intense gets removed first)
      remove_order = length(which(pks$peakgr == peakgr)):1
    )
    dt$P <- p_vals[match(dt$peakid, pks$peakid)]
    dt
  })
  peaks <- do.call("rbind", peaks)

  ## save probabilities
  write.csv(peaks, paste0(out_dir, "probs_missingness.csv"), row.names = FALSE)
}


# helpers for modelMISS -------------------------------------------------------------------------------------------
## helper functions for modelMISS function

## optimiseBETA function finds intercept beta0 by numerically solving the equation defined in Do et al. 2018
optimiseBETA <- function(n, b1, miss, ints) {
  out <- data.frame()
  b0 <- -100
  p_sum <- -Inf
  while (p_sum < miss) {
    p_vals <- sapply(1:n, FUN = getPROB, ints = ints, b0 = b0, b1 = b1)
    if (!NaN %in% p_vals) {
      p_sum <- sum(p_vals)
      out <- rbind(out, data.frame(b0 = b0, p_sum = p_sum))
    } else {
      p_sum <- -Inf
    }
    b0 <- b0 + 1
  }
  ## get more precise b0
  b0 <- b0 - 1
  out <- data.frame()
  p_sum <- -Inf
  while (p_sum < miss) {
    p_vals <- sapply(1:n, FUN = getPROB, ints = ints, b0 = b0, b1 = b1)
    if (!NaN %in% p_vals) {
      p_sum <- sum(p_vals)
      out <- rbind(out, data.frame(b0 = b0, p_sum = p_sum))
    } else {
      p_sum <- -Inf
    }
    b0 <- b0 + 0.001
  }
  b0 <- out[ifelse(nrow(out) > 1, nrow(out) - 1, 1), "b0"]
  return(b0)
}

## function
getPROB <- function(ints, p, b0, b1) {
  i <- ints[p]
  a <- b0 + b1 * i
  logistic(a)
}

## logistic function of the coefficient and intercept values
logistic <- function(a) {
  exp(a) / (1 + exp(a))
}

## scale peaks intensities to enable logistic function calculation
scalePEAKS <- function(ints, from = 1e-20, to = 1) {
  (ints - min(ints)) / max(ints - min(ints)) * (to - from) + from
}
