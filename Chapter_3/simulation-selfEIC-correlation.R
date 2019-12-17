## -------------------------------------------------------------------------------------------
## Fuctions used for a simulation of EIC correlation between identical peaks
## -------------------------------------------------------------------------------------------

# checkPk -----------------------------------------------------------------
## Helper function to evaluate how well a centWave detected peak fits Gaussian shape
## Function fits Nonlinear Least Squares to the chromatographic profile of the centWave detected peak using xcms-derived Gaussian parameters
checkPk <- function(pk, # index of the single peak in the EIC list
                    eic # EIC list from which a single peak is to be extracted from
) {
  eic_pk <- data.frame(
    scan.raw = as.numeric(gsub("F1.S", "", names(MSnbase::rtime(eic[[pk]])))),
    rt = MSnbase::rtime(eic[[pk]]),
    into = MSnbase::intensity(eic[[pk]])
  )
  pos <- min(eic_pk$scan.raw):max(eic_pk$scan.raw)
  eic_pk$scan <- match(eic_pk$scan.raw, pos)

  ## use xcms-derived Gaussian parameters for the peak
  params <- c("mu" = pks[pk, "mu"], "sigma" = pks[pk, "sigma"], "h" = pks[pk, "h"])
  result <- tryCatch(
    {
      nls(into ~ k * exp(-1 / 2 * (scan.raw - mu)^2 / sigma^2), start = c(mu = params[[1]], sigma = params[[2]], k = params[[3]]), data = eic_pk)
    },
    error = function(er) {
      return(er)
    }
  )
  if (class(result)[1] != "nls") {
    ## if peak is so far away from the gaussian shape that NLS does not get fit
    return(NULL)
  }
  fit <- result

  ## update parameters
  params <- summary(fit)$parameters[, "Estimate"]

  plot(into ~ scan.raw,
    data = eic_pk, type = "b", cex = .5,
    xaxt = "n", xlab = "Scan number",
    ylab = "Intensity",
    main = paste0("Peak no: ", pk)
  )
  plot(function(x) params[3] * exp(-1 / 2 * (x - params[1])^2 / params[2]^2),
    add = T,
    xlim = range(eic_pk$scan.raw),
    col = "red"
  )
  axis(1, at = eic_pk$scan.raw, labels = eic_pk$scan)
  return(fit)
}

# corPk -------------------------------------------------------------------
## take a characterically-looking peak from a representative sample
## obtain Gaussian parameters (as provided by centWave peak-picker with argument fitgauss = TRUE)
## observe self-EIC-correlation as apexes move further away
corPk <- function(pk, # index of the single peak in the EIC list
                  eic, # EIC list from which a single peak is to be extracted from
                  out_dir # path to output directory
) {
  require(ggplot2)
  ## extract chromatogram from selected peak
  tab <- data.frame(
    scan.raw = as.numeric(gsub("F1.S", "", names(MSnbase::rtime(eic[[pk]])))),
    rt = MSnbase::rtime(eic[[pk]]),
    into = MSnbase::intensity(eic[[pk]])
  )
  tab$scan <- sapply(tab$scan.raw, function(x) {
    match(x, min(tab$scan.raw):max(tab$scan.raw))
  })

  ## use xcms-derived Gaussian parameters for the peak
  params <- c("mu" = pks[pk, "mu"], "sigma" = pks[pk, "sigma"], "h" = pks[pk, "h"])
  tab_params <- data.frame(
    value = c(
      params["mu"],
      params["mu"] + params["sigma"],
      params["mu"] - params["sigma"],
      params["mu"] + (2 * params["sigma"]),
      params["mu"] - (2 * params["sigma"])
    ),
    param = c(
      "mu",
      rep("sigma", 4)
    ),
    stringsAsFactors = F,
    row.names = NULL
  )

  #### ---- (A) plot chrom peak and fitted Gaussian
  cols <- viridis::viridis(
    n = 3,
    end = 0.8, # omitting bright yellow
    alpha = 0.8
  )
  g <- ggplot(
    data = tab,
    aes(x = scan.raw, y = into)
  ) +
    geom_point() +
    geom_vline(
      data = tab_params,
      aes(xintercept = value, col = param)
    ) +
    stat_function(
      fun = function(x) params[3] * exp(-1 / 2 * (x - params[1])^2 / params[2]^2),
      geom = "line", color = "grey"
    ) +
    scale_color_manual(
      values = cols,
      labels = c("Median", "+/- 1/2 SD", "Fitted Gaussian"),
      name = ""
    ) +
    scale_x_continuous(
      name = "Scan",
      breaks = tab$scan.raw[seq(0, nrow(tab), by = 5)],
      labels = tab$scan[seq(0, nrow(tab), by = 5)]
    ) +
    ylab("Intensity") +
    theme_bw() +
    theme(legend.position = "bottom")

  grid::grid.newpage()
  pdf(
    width = 6, height = 6,
    file = paste0(out_dir, "/", sname, "_checkPeak_pk", pk, "_chrom", ".pdf")
  )
  grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
  dev.off()

  #### ---- (B) get self-correlation
  scans <- tab$scan
  # find where apex of the peak is
  sc_apex <- scans[which(abs(params["mu"] - tab$scan.raw) == min(abs(params["mu"] - tab$scan.raw)))]
  # find maximum scan no after moving
  sc_max <- nrow(tab) * 2 + 1
  sc_seq <- c(0, scans)
  sc_min <- 1
  out <- data.frame()
  for (sc in sc_seq) {
    # move apex by the sc number of scans
    tab$scan_now <- tab$scan + sc
    common_scan <- base::intersect(tab$scan[which(!is.na(tab$into))], tab$scan_now[which(!is.na(tab$into))])
    if (length(common_scan) > 3) {
      x <- tab$into[which(tab$scan %in% common_scan)]
      y <- tab$into[which(tab$scan_now %in% common_scan)]
      cc <- cor(x, y, method = "pearson")
      g <- ggplot(tab) +
        # area of scans in common
        geom_rect(aes(xmin = min(common_scan) - 0.5, xmax = max(common_scan) + 0.5, ymin = -Inf, ymax = Inf), alpha = 0.01, fill = "grey") +
        # original peak
        geom_segment(aes(x = scan, xend = scan, y = 0, yend = into), color = cols[1], alpha = 0.6) +
        geom_point(aes(y = into, x = scan), color = cols[1], alpha = 0.6) +
        # moved peak
        geom_point(aes(y = into, x = scan_now), color = cols[2], alpha = 0.6) +
        geom_segment(aes(x = scan_now, xend = scan_now, y = 0, yend = into), color = cols[2], alpha = 0.6) +
        scale_x_continuous(limits = c(0, sc_max + 0.5)) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(colour = "black", size = 0.5)
        ) +
        xlab("Scan") +
        ylab("Intensity") +
        ggtitle(paste0("Scan window: ", sc, ". Correlation: ", round(cc, digits = 3)))

      grid::grid.newpage()
      pdf(
        width = 6, height = 6,
        file = paste0(out_dir, "/", sname, "_checkPeak_pk", pk, "_scanWindow-", sc, ".pdf")
      )
      grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
      dev.off()
    } else {
      message(paste("Peaks do not overlap enough. Scan:", sc))
      cc <- NA
    }
    out_tab <- data.frame(scan = sc, cor = cc)
    out <- rbind(out, out_tab)
  }

  #### ---- (C) plot cor coeficient vs distance (scans and sd)
  out <- rbind(out, out)
  # add with -minus scan cor to make symmetrical graph
  out$scan[(length(sc_seq) + 1):nrow(out)] <- -out$scan[(length(sc_seq) + 1):nrow(out)]

  g <- ggplot(data = out, aes(x = scan, y = cor)) +
    geom_point(color = "grey", na.rm = T) +
    geom_point(
      data = out %>% filter(cor > 0.9),
      aes(x = scan, y = cor, color = "Correlation > 0.9")
    ) +
    scale_color_manual(values = cols[1], name = "", labels = "Correlation > 0.9") +
    xlab("Apex distance, scans") +
    ylab("Correlation") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(colour = "black", size = 0.5)
    ) +
    theme(legend.position = "bottom")

  grid::grid.newpage()
  pdf(
    width = 6, height = 6,
    file = paste0(out_dir, "/", sname, "_checkPeak_pk", pk, "_correlation", ".pdf")
  )
  grid::grid.draw(rbind(ggplotGrob(g), size = "last"))
  dev.off()
}
