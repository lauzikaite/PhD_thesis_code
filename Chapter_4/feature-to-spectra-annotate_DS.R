## -------------------------------------------------------------------------------------------
## Fuctions used for feature-to-spectra matching
## -------------------------------------------------------------------------------------------

# annotate_DS -------------------------------------------------------------
## main function facilitating dataset annotation with a database table
annotate_DS <- function(matrix_file, # path to features table
                        db, # database table
                        mz_err, # mz error
                        rt_err, # rt error
                        major_spec = TRUE, # should only major rt specs be used in annotation
                        thrs, # confidence score threshold
                        save_plot = FALSE, # should plots be generated
                        out_dir # path to output directory
) {
  require(ggplot2)
  require(gridExtra)
  options(warn = -1)

  if (file.access(names = out_dir, mode = 2) != 0) {
    stop(
      "You don't have permission to write in the directory: ", out_dir,
      "\n Select out_dir with read-write permissions"
    )
  }

  # set default mz/rt window for matching
  if (missing(mz_err)) {
    mz_err <- 0.01
    message("mz_err value not selected, using default: 0.01")
  }

  if (missing(rt_err)) {
    rt_err <- 0.2
    message("rt_err value not selected, using default: 0.2")
  }

  # set default MatschScore threshold
  if (missing(thrs)) {
    thrs <- 0
    message("thrs value not selected, using default: 0")
  }

  # load features
  ds <- read.csv(matrix_file, header = TRUE, stringsAsFactors = FALSE)
  req_cnames <- c("m.z", "Retention.Time")
  if (any(!req_cnames %in% colnames(ds))) {
    stop("matrix_file must contain columns 'm.z' and 'Retention.Time'!")
  }
  # add features ID column
  if (!"ID" %in% colnames(ds)) {
    ds$ID <- seq(nrow(ds))
  }

  # make a copy of database
  db_anno <- db

  ## if only major rt groups should be considered in matching
  if (major_spec == TRUE) {
    message("only major rt groups in the chemical standards will be considered")
    db_anno <- db_anno[db_anno$major.spec == TRUE, ]
  }

  # calculate mz and rt for all spectra in the databse
  db_anno$MZ_high <- db_anno[, "MZ"] + mz_err
  db_anno$MZ_low <- db_anno[, "MZ"] - mz_err

  db_anno$RT_high <- db_anno[, "spec.rt"] + rt_err
  db_anno$RT_low <- db_anno[, "spec.rt"] - rt_err

  out <- rep(list(list()), length(unique(db_anno$chem.NPCID)))
  names(out) <- as.character(unique(db_anno$chem.NPCID))
  out_thrs <- data.frame()

  # get progress bar
  pb <- txtProgressBar(min = 0, max = length(unique(db_anno$chem.NPCID)), style = 3)

  ##### for each DB compound, extract rt groups ----
  for (i in unique(db_anno$chem.NPCID)) {

    # extract compound dataframe
    compound <- subset(db_anno, chem.NPCID == i)
    compound$ds <- rep(FALSE, nrow(compound)) # redundant, required for ggplot df build

    # extract rt groups
    rtg <- unique(compound$spec.no)

    # create three df for each compound: (1) for original db; (2) for ds matches (if any); (3) for int scores for DB and matches (if any)
    out[[i]] <- list(
      "db" = data.frame(compound),
      "ds" = data.frame(),
      "sc" = data.frame()
    )

    ##### for each rt group, find matches ----
    mat_cen <- lapply(rtg, function(r) {
      rtgroup <- subset(compound, spec.no == r)
      peaks <- nrow(rtgroup)
      mat <- data.frame()

      # to debug issue associated with db built with higher threshold, which can remove rtgroups
      if (peaks > 0) {

        ##### for each peak in rt group, find matches ----
        for (p in (unique(rtgroup$peak.no))) {
          peakgroup <- subset(rtgroup, peak.no == p)

          mz_h <- ds$"m.z" <= rtgroup[p, "MZ_high"]
          mz_l <- ds$"m.z" >= rtgroup[p, "MZ_low"]

          rt_h <- ds$"Retention.Time" <= rtgroup[p, "RT_high"]
          rt_l <- ds$"Retention.Time" >= rtgroup[p, "RT_low"]

          sc <- mz_h + mz_l + rt_h + rt_l
          matches <- which(sc == 4)

          # does this peak have any matches in the dataset?
          if (length(matches) > 0) {
            matches <- data.frame(
              "ID" = ds[matches, "ID"],
              "MZ" = ds[matches, "m.z"],
              "spec.rt" = ds[matches, "Retention.Time"],
              "spec.norm.it" = rep(peakgroup$spec.norm.it, length(matches)),
              "chem.norm.it" = rep(peakgroup$chem.norm.it, length(matches)),
              "spec.no" = rep(r, length(matches)),
              "major.spec" = rep(peakgroup$major.spec, length(matches)),
              "peak.no" = rep(p, length(matches)),
              "chem.NPCID" = rep(peakgroup$chem.NPCID, length(matches)),
              "chem.name" = rep(peakgroup$chem.name, length(matches)),
              "ds" = rep(TRUE, length(matches)), # redundant, required for ggplot df build
              row.names = NULL
            )
            mat <- rbind(mat, matches)
          }
        }


        ##### if rt group has A MATCH in the dataset ----
        if (nrow(mat) > 0) {

          #### calculate distance to central rt for each match ----
          # (1) if >2 matches that are NOT the same rt, find DS rt density:
          if (nrow(mat) > 2 & !all(mat$spec.rt == mat$spec.rt[1])) {
            # density weighted by DB peak intensity
            it_weight <- mat$chem.norm.it / sum(mat$chem.norm.it)
            d <- density(mat$spec.rt, weights = it_weight)
            c <- d$x[which(d$y == max(d$y))]
            if (length(c) > 1) {
              c <- median(c)
            }
            cen <- data.frame(
              "den" = rep(TRUE, nrow(mat)),
              "cen" = rep(c, nrow(mat)),
              check.names = FALSE
            )
            mat <- cbind(mat, cen)
          } else {

            # (2) if >2 matches, but they are the SAME, use DS rt:
            if (nrow(mat) >= 2 & all(mat$spec.rt == mat$spec.rt[1])) {
              c <- mat$spec.rt[1]
              cen <- data.frame(
                "den" = rep(FALSE, nrow(mat)),
                "cen" = rep(c, nrow(mat)),
                check.names = FALSE
              )
              mat <- cbind(mat, cen)
            } else {
              # (3) if <2 matches, use original DS rt:
              c <- unique(mat$spec.rt)
              cen <- data.frame(
                "den" = rep(FALSE, nrow(mat)),
                "cen" = rep(c, nrow(mat)),
                check.names = FALSE
              )
              mat <- cbind(mat, cen)
            }
          }
        }
      }
    })
    mat <- do.call(rbind, mat_cen)

    #### if COMPOUND has a match, CONTINUE ----
    # update progress bar
    y <- which(unique(db_anno$chem.NPCID) == i)
    setTxtProgressBar(pb, y)

    if (is.null(nrow(mat))) {
      next
    }

    mat$d <- mat$spec.rt - mat$cen

    #### find the TOP match for each peak, rt group-wise ----
    mat_top <- lapply(unique(mat$spec.no), function(g) {
      tmp <- subset(mat, spec.no == g)

      # calc standard deviation
      if (nrow(tmp) == 1) {
        std <- 0
      } else {
        std <- sd(tmp$spec.rt)
      }
      # RT window is +/- 2 standard deviations from the RT center
      std_l <- unique(tmp[, "cen"]) - abs(std) * 2
      std_h <- unique(tmp[, "cen"]) + abs(std) * 2

      tmp[, "cen_l"] <- std_l
      tmp[, "cen_h"] <- std_h

      # for each peak in the rt groug
      mat_top <- lapply(unique(tmp$peak.no), function(p) {
        tmp <- subset(tmp, peak.no == p)

        # if a single match, check if is within [2*std] from the cen
        if (nrow(tmp) == 1) {
          if (with(tmp, spec.rt <= std_h & spec.rt >= std_l)) {
            tmp$top <- TRUE
          } else {
            tmp$top <- FALSE
          }
          mat_top <- tmp
        } else {
          # if more than one match, find the closest match
          tmp$top <- rep(FALSE, nrow(tmp))
          t <- which(order(abs(tmp$d), decreasing = F) == 1)
          tmp[t, "top"] <- TRUE
          mat_top <- tmp
        }
      })

      mat_top <- do.call(rbind, mat_top)
    })

    mat <- do.call(rbind, mat_top)
    out[[i]]$ds <- mat

    #### if the MAJOR rtg group have a match, calc MatchScore ----
    if (any(out[[i]]$ds$major.spec == TRUE)) {

      #### calc  matching score for TOP matches ----
      db_int <- sum(subset(out[[i]]$db, major.spec == TRUE)$chem.norm.it)
      ds_int <- sum(subset(out[[i]]$ds, major.spec == TRUE & top == TRUE)$chem.norm.it)
      db_peaks <- length(subset(out[[i]]$db, major.spec == TRUE)$peak.no)
      ds_peaks <- length(subset(out[[i]]$ds, major.spec == TRUE & top == TRUE)$peak.no)
      sc <- data.frame(
        "DBInt" = db_int,
        "DSInt" = ds_int,
        "DSIntSc" = (ds_int / db_int),
        "DSPeakIntSc" = (ds_int / db_int) * ds_peaks,
        "DBMax" = 1 * db_peaks,
        "MatchScore" = round(((ds_int / db_int) * ds_peaks) / (1 * db_peaks), digits = 3)
      )
      out[[i]]$sc <- sc

      #### if MatchScore ABOVE thres, CONTINUE ----
      if (nrow(out[[i]]$sc) > 0) { # redundant, to account for an empty sc dataframe if score is NA
        if (out[[i]]$sc$MatchScore > thrs) {

          #### Write to output ----
          mat_thrs <- subset(out[[i]]$ds, top == TRUE)[, c(
            "ID",
            "MZ",
            "spec.rt",
            "spec.norm.it",
            "chem.norm.it",
            "spec.no",
            "major.spec",
            "peak.no",
            "chem.NPCID",
            "chem.name"
          )]
          mat_thrs <- cbind(mat_thrs, out[[i]]$sc)
          out_thrs <- rbind(out_thrs, mat_thrs)

          # for each of the selected match, extract corresponding compound's peak
          out[[i]]$db$top <- rep(F, nrow(out[[i]]$db))
          sel <- which(with(out[[i]]$db, spec.no %in% mat_thrs$spec.no & peak.no %in% mat_thrs$peak.no))
          out[[i]]$db[sel, "top"] <- rep(T, length(sel))

          #### PLOTS: prep ----
          if (save_plot == TRUE) {
            plot_dat <- merge(out[[i]]$db,
              out[[i]]$ds,
              by = intersect(
                colnames(out[[i]]$db),
                colnames(out[[i]]$ds)
              ),
              all = TRUE
            )
            plot_rtg <- lapply(unique(subset(plot_dat, ds == FALSE)$spec.no), function(n) {
              unique(subset(plot_dat, ds == FALSE & spec.no == n, select = c("RT_low", "RT_high")))
            })
            plot_rtg <- do.call(rbind, plot_rtg)
            rt_min <- min(plot_rtg$RT_low)
            rt_max <- max(plot_rtg$RT_high)
            plot_cen <- lapply(unique(subset(plot_dat, ds == TRUE)$spec.no), function(n) {
              unique(subset(plot_dat, ds == TRUE & spec.no == n, select = c("cen_l", "cen_h")))
            })
            plot_cen <- do.call(rbind, plot_cen)

            #### PLOTS(A) plot DB vs DS ----
            p_regions <- ggplot() +
              # first add areas for rt windows and stdev to leave them in the first layer
              geom_rect(data = plot_rtg, aes(xmin = RT_low, xmax = RT_high, ymin = -Inf, ymax = Inf), alpha = 0.08, fill = "grey") +
              geom_rect(data = plot_cen, aes(xmin = cen_l, xmax = cen_h, ymin = -Inf, ymax = Inf), alpha = 0.1, fill = "#b2df8a") +
              theme_bw(base_size = 14) +
              scale_x_continuous(limits = c(rt_min, rt_max), name = "Retention Time") +
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "none",
                axis.title.x = element_blank(),
                axis.title.y = element_blank()
              ) +
              geom_vline(data = subset(plot_dat, den == TRUE), aes(xintercept = cen), size = 0.2, colour = "black", linetype = "longdash")

            p1 <-
              p_regions +
              ggtitle(
                paste0(
                  unique(out[[i]]$db$chem.NPCID),
                  ", confidence score: ", sc$MatchScore
                )
              ) +
              # plot OTHER rt groups: DB only
              geom_point(data = subset(plot_dat, major.spec == FALSE), aes(spec.rt, MZ, size = chem.norm.it, shape = ds), alpha = 0.6, colour = "black") +
              # plot MAJOR rt: shape by ds/db, size by int, color by top
              geom_point(data = subset(plot_dat, major.spec == TRUE), aes(spec.rt, MZ, size = chem.norm.it, colour = as.character(top), shape = ds), alpha = 0.6) +
              scale_shape_manual(
                name = "",
                labels = c("DB", "Matches"),
                values = c(19, 8)
              ) +
              scale_size_continuous(
                guide = FALSE,
                range = c(2, 6)
              ) +
              scale_colour_manual(
                name = "",
                labels = c("Not selected", "Selected match"),
                values = setNames(c("#35B779FF", "#440154FF"),
                  nm = c(TRUE, FALSE)
                )
              ) +
              theme(
                plot.title = element_text(hjust = 0.5),
                legend.position = "top"
              )

            #### PLOTS: (B) if DEN estimate for the MAJOR rt group is available, make second plot ----
            if (any(subset(plot_dat, major.spec == TRUE & ds == TRUE)$den == TRUE)) {
              plot_dat <- subset(plot_dat, ds == TRUE & den == TRUE & major.spec == TRUE)
              w <- plot_dat$chem.norm.it / sum(plot_dat$chem.norm.it)
              p2 <-
                p_regions +
                geom_density(
                  data = plot_dat, aes(spec.rt, y = ..scaled.., group = spec.no, fill = as.character(major.spec), weight = w),
                  alpha = 0.5, size = 0.01, fill = "#35B779FF"
                )
              p1_g <- ggplotGrob(p1)
              p2_g <- ggplotGrob(p2)
              grid::grid.newpage()
              pp <- gridExtra:::gtable_rbind(p1_g, p2_g, size = "max")
            } else {
              pp <- p1
            }
            suppressMessages(ggsave(
              filename = file.path(out_dir, paste0(i, ".png")),
              plot = pp,
              height = 200, width = 200, units = "mm",
              dpi = 200
            ))
          }
        }
      }
    }
  }

  #### FINAL output ----
  # return concatinated list of annotated features
  features <- out_thrs[, match(
    c("ID", "MZ", "spec.rt", "chem.NPCID", "chem.name", "MatchScore"),
    colnames(out_thrs)
  )]
  # extract annotations for every feature
  features_anno <- lapply(unique(features$ID),
    FUN = extract_ANNO,
    features = features,
    id_name = "ID"
  )
  features_anno <- do.call("rbind", features_anno)
  features_final <- merge(ds,
    features_anno,
    by = "ID",
    all = TRUE
  )
  # extract annotations for every chemical compound
  chems_anno <- lapply(unique(db$chem.NPCID),
    FUN = extract_ANNO,
    features = features,
    id_name = "chem.NPCID"
  )
  chems_anno <- do.call("rbind", chems_anno)
  chems_anno <- chems_anno[chems_anno$anno.count > 0, ] # retain only annotated compounds
  chems_anno$MatchScore <- sapply(out[chems_anno$chem.NPCID], function(x) {
    x[["sc"]]$MatchScore
  })
  chems_final <- merge(unique(db[, c("chem.NPCID", "chem.name")]),
    chems_anno,
    by = "chem.NPCID"
  )
  # write output
  write.csv(features_final,
    file = file.path(out_dir, "DS_annotated.csv"),
    row.names = FALSE,
    col.names = TRUE
  )
  write.csv(chems_final,
    file = file.path(out_dir, "DB_annotated.csv"),
    row.names = FALSE,
    col.names = TRUE
  )
}


# extract_ANNO ------------------------------------------------------------
## helper function to extract annotations for either peaks in the dataset or chemicals in the database
extract_ANNO <- function(id_value, # iterator value, either chem.NPCID or dataset peak ID
                         features, # annotation results table
                         id_name # character equal to either chem.NPCID or ID
) {
  req_cnames <-
    if (id_name == "chem.NPCID") {
      c("ID")
    } else {
      c("chem.NPCID", "chem.name", "MatchScore")
    }
  dat <- features[features[, match(id_name, colnames(features))] == id_value, ]
  dat <- dat[order(dat$MatchScore, decreasing = TRUE), ] # order by confidence score
  dat_by_id <- unique(as.data.frame(dat[, req_cnames]))
  cbind(
    setNames(data.frame(id_value,
      nrow(dat_by_id),
      stringsAsFactors = FALSE
    ),
    nm = c(id_name, "anno.count")
    ),
    setNames(as.data.frame(lapply(seq(req_cnames), function(ncol) {
      paste(dat_by_id[, ncol], collapse = "|")
    })), nm = req_cnames)
  )
}
