## -------------------------------------------------------------------------------------------
## Fuctions used for spectral scaling analysis
## -------------------------------------------------------------------------------------------
## the following packages must be installed
## these packages do not need to be loaded
# require(cowplot)
# require(viridis)
# require(dplyr)
# require(magrittr) # for pipe only
# require(grid)
# require(gridExtra)
# require(tidyr)

# getCOS ------------------------------------------------------------------
## obtain cosine betweem two PCS: target and matches peaks
getCOS <- function(spec, target_peaks, matched_peaks, scale, norm) {
  ## build vector for target peaks
  target_vec <- getVEC(spec = spec, peaks = target_peaks, scale = scale, norm = norm)
  ## build vector for matches peaks
  matched_vec <- getVEC(spec = spec, peaks = matched_peaks, scale = scale, norm = norm)
  ## calculate cosine of the angle between the vectors
  cos_norm <- (sum(target_vec * matched_vec)) / ((sqrt(sum(
    target_vec * target_vec
  ))) * (sqrt(sum(
    matched_vec * matched_vec
  ))))
  return(cos_norm)
}

# getVEC ------------------------------------------------------------------
## build vector using provided peaks m/z and intensity values
## using reference spec to find where the peaks m/z values belong to
getVEC <- function(spec, peaks, scale, norm) {
  peaks$bin <- findInterval(x = peaks$mz, spec$mz)
  peaks <- peaks[which(!duplicated(peaks$bin)), c("mz", "bin", "into")]
  spec$into <- 0
  
  if (scale == "weight") {
    ## scaling by the weight of the features, this is not going to be normalised
    peaks$into <- apply(peaks[, c("mz", "into")],
                        1,
                        FUN = scaleWEIGHT
    )
  } else {
    if (scale == "sqrt") {
      ## scaling by taking the sqrt
      peaks$into <- sqrt(peaks$into)
    }
  }
  spec[match(peaks$bin, spec$bin), "into"] <- peaks$into
  
  if (norm == TRUE) {
    ## normalising by the magnitude of the spectral vector
    vector <- normMAGN(into = spec$into)
  } else {
    ## no normalisation
    vector <- spec$into
  }
  return(vector)
}

scaleWEIGHT <- function(x) {
  x[["into"]]^0.6 * x[["mz"]]^3
}

# normMAGN ----------------------------------------------------------------
## normalise scaled vector to the total magnitude of the spectral vector
normMAGN <- function(into) {
  into / (sqrt(sum(into * into)))
}

# plotSCENARIO ------------------------------------------------------------
## Plot target and matched PCS as vectors
## Make three sub-plots: one for each scaling method
plotSCENARIO <- function(scen, peakgroups_dat, scen_name, ref, id) {
  ## select a sample to compare with the reference spec
  run_order_sel <- ifelse(nrow(scen) == 1, scen$run_order, sample(x = c(scen$run_order), size = 1))
  ds_spec <- peakgroups_dat[peakgroups_dat$run_order == run_order_sel, ]
  
  ## prep objects for plotting
  spec_cos <- subset(scen, run_order == run_order_sel)
  gg_cols <- c("black", viridis::viridis(length(unique(ref$adduct)) + 1, begin = 0.1)[2:(length(unique(ref$adduct)) + 1)])
  gg_labels <- setNames(c(
    paste0("Unscaled: cos ", round(spec_cos$no_scale_norm, digits = 3)),
    paste0("Sqrt-scaled: cos ", round(spec_cos$scale_sqrt_norm, digits = 3)),
    paste0("Weight-scaled: cos ", round(spec_cos$scale_weight_norm, digits = 3))
  ),
  nm = c("none", "sqrt", "weight")
  )
  
  gg_unscaled <- plotSPECTRA(
    dt = rbind(
      ref %>%
        mutate(spec = "ref"),
      ds_spec %>%
        mutate(spec = "ds")
    ) %>%
      mutate(scaling = "none") %>%
      dplyr::select(mz, rt, run_order, adduct, scaling, into = into, spec),
    gg_cols = gg_cols,
    gg_labels = gg_labels
  )
  gg_sqrt <- plotSPECTRA(
    dt = rbind(
      ref %>%
        mutate(spec = "ref"),
      ds_spec %>%
        mutate(spec = "ds")
    ) %>%
      mutate(scaling = "sqrt") %>%
      dplyr::select(mz, rt, run_order, adduct, scaling, into = into_sqrt_norm, spec),
    gg_cols = gg_cols,
    gg_labels = gg_labels
  )
  gg_weight <- plotSPECTRA(
    dt = rbind(
      ref %>%
        mutate(spec = "ref"),
      ds_spec %>%
        mutate(spec = "ds")
    ) %>%
      mutate(scaling = "weight") %>%
      dplyr::select(mz, rt, run_order, adduct, scaling, into = into_weight_norm, spec),
    gg_cols = gg_cols,
    gg_labels = gg_labels
  )
  
  ## to make plots with equal widths (y-axis labels have different widths), rbind them with size = "first"
  gg_grobs <- lapply(
    list(
      gg_unscaled +
        theme(legend.position = "none"),
      gg_sqrt +
        theme(legend.position = "none"),
      gg_weight +
        theme(legend.position = "none")
    ),
    ggplotGrob
  )
  gg_rbind <- do.call(rbind, c(gg_grobs, size = "first"))
  gg_legend <- cowplot::get_legend(gg_weight +
                                     theme(legend.position = "bottom"))
  gg_title <- cowplot::get_title(gg_unscaled +
                                   ggtitle(paste0(standards$cpdName[id_inds][1], ", ", id)) +
                                   theme(plot.title = element_text(hjust = 0.5)))
  
  ## add title and legend
  grid::grid.newpage()
  gg <- gridExtra::arrangeGrob(
    gg_title,
    gg_rbind,
    gg_legend,
    nrow = 3,
    heights = c(1, 10, 1)
  )
  ggsave(paste0("spectralSimilarity_scen", scen_name, "_", id, ".png"),
         plot = gg, path = out_dir, height = 250, width = 150, units = "mm", dpi = "print"
  )
}

# plotSPECTRA -------------------------------------------------------------
## helper function for plotSCEN
plotSPECTRA <- function(dt, gg_cols, gg_labels) {
  ggplot(data = dt) +
    ## make top spectra: 1st spectra in the pair
    geom_segment(
      data = subset(dt, spec == "ref"),
      aes(x = mz, xend = mz, y = 0, yend = into, color = as.factor(adduct)),
      size = 1, na.rm = TRUE
    ) +
    ## make bottom spectra: 2nd spectra in the pair
    geom_segment(
      data = subset(dt, spec == "ds"),
      aes(x = mz, xend = mz, y = 0, yend = -into, color = as.factor(adduct)),
      size = 1, na.rm = TRUE
    ) +
    scale_color_manual(
      name = "Adducts",
      values = gg_cols,
      labels = c("Un-assigned", standards[id_inds, "ROI_mz"])
    ) +
    scale_x_continuous("m/z") +
    scale_y_continuous("Scaled intensity", limits = c(-max(dt$into), max(dt$into))) +
    facet_wrap(~scaling, ncol = 1, scales = "free", labeller = as_labeller(gg_labels)) +
    geom_hline(data = data.frame(y = 0), aes(yintercept = y)) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.1)
    )
}
