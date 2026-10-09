## =============================================================================
## nuclei_sybr_peak_analysis.R
##
## SYBR Green I nuclei cytometry: FACS density panels, peak fitting, bootstrap intervals.
##
## Steps (run top to bottom):
##   1. read the FCS files
##   2. FACS panels (FSC vs SYBR, density-coloured, SYBR+ gate) for all selected files
##   3. for each stained file: find the peaks, fit each one in a window (Gaussian on log10
##      + baseline), measure the fold change between neighbouring peaks
##   4. bootstrap the peak positions and fold changes
##   5. save figures (PDF) and tables (CSV) to results/
##
## Usage (from the project folder):   source("nuclei_sybr_peak_analysis.R")
##                                or: Rscript nuclei_sybr_peak_analysis.R
## Requires: flowCore (Bioconductor), ggplot2, patchwork, scales. Optional: ggrastr.
## =============================================================================

## ---- 0. Packages ------------------------------------------------------------
pkgs <- c("flowCore", "ggplot2", "patchwork", "scales")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("Missing packages: ", paste(miss, collapse = ", "),
       "\n  install.packages(c(\"ggplot2\", \"patchwork\", \"scales\"))",
       "\n  BiocManager::install(\"flowCore\")")
suppressPackageStartupMessages({ library(flowCore); library(ggplot2); library(patchwork) })

## ---- 1. Configuration -------------------------------------------------------
fcs_dir           <- "raw_fcs"                       # searched recursively
file_pattern      <- "tmuris-female.*\\.fcs$"        # which files to read (regex on file name)
unstained_pattern <- "unstained"                     # files matching this are controls
sybr_ch           <- "FL53-A"                        # SYBR Green I, 507/19 detector
fsc_ch            <- "FSC04-A"                       # 488 nm forward scatter
out_dir           <- "results"

gate      <- 10            # SYBR+ gate; events above it are analysed
x_range   <- c(1, 1e5)     # SYBR axis, log10: 10^0 - 10^5
ylab_max  <- 512           # relabel FSC ticks to a 0-512 display scale (labels only); NULL = raw values
rasterise_points <- FALSE  # TRUE: dots as an image inside the PDF (needs ggrastr), smaller files

# peak detection (apexes of a smoothed density of events above peak_floor)
peak_floor <- 20; peak_adjust <- 0.5; peak_min_rel <- 0.1
# optional: fix the peak positions for a file instead of auto-detecting, e.g.
# centers_override <- list("tmuris-female-stained_TMURF v2 Run 1 20231101162939.fcs" = c(27.9, 72.8, 140.7))
centers_override <- list()

n_boot <- 500              # bootstrap resamples
seed   <- 1

red <- "#C0392B"
pal_ref <- colorRampPalette(c("#6A7BE0", "#2E3FD0", "#00B7EB", "#2EBF5C",
                              "#F2E600", "#FF8A00", "#E8000B"))(100)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 2. Helper functions ----------------------------------------------------

# Save a figure. PDFs use the Cairo device so Unicode symbols and fonts come out correctly.
save_fig <- function(filename, plot, width, height, ...) {
  if (grepl("\\.pdf$", filename, ignore.case = TRUE)) {
    dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
    ggsave(filename, plot, width = width, height = height, device = dev, ...)
  } else {
    ggsave(filename, plot, width = width, height = height, dpi = 300, ...)
  }
}

read_events <- function(file) {
  ff <- read.FCS(file, transformation = FALSE, truncate_max_range = FALSE)
  m  <- exprs(ff); cn <- unname(colnames(m))
  if (!all(c(sybr_ch, fsc_ch) %in% cn))
    stop("channels ", sybr_ch, "/", fsc_ch, " not found in ", basename(file),
         ". Available: ", paste(cn, collapse = ", "))
  d <- data.frame(sybr = as.numeric(m[, match(sybr_ch, cn)]),
                  fsc  = as.numeric(m[, match(fsc_ch,  cn)]))
  list(d = d[is.finite(d$sybr) & is.finite(d$fsc), ], file = file)
}

# Events per grid cell (nb x nb over the plotted area): used to colour dots by local density
point_density <- function(d, xlim, ymax, nb = 200, what = c("event", "cell")) {
  what <- match.arg(what)
  tx <- pmin(pmax(log10(pmax(d$sybr, 1e-12)), log10(xlim[1])), log10(xlim[2]))
  ty <- pmin(pmax(d$fsc, 0), ymax)
  ix <- findInterval(tx, seq(log10(xlim[1]), log10(xlim[2]), length.out = nb + 1), all.inside = TRUE)
  iy <- findInterval(ty, seq(0, ymax, length.out = nb + 1), all.inside = TRUE)
  cell <- (iy - 1L) * nb + ix
  cnt  <- tabulate(cell, nbins = nb * nb)
  if (what == "cell") cnt[cnt > 0] else cnt[cell]
}

# FACS panel: FSC vs SYBR, density-coloured dots, green SYBR+ gate with the percentage
plot_facs_panel <- function(e, title = NULL, ymax, cmax, nb = 200, pt_size = 0.25) {
  d   <- e$d
  pct <- 100 * mean(d$sybr > gate)
  d$dens <- point_density(d, x_range, ymax, nb)
  d$sybr <- pmin(pmax(d$sybr, x_range[1]), x_range[2])   # log axis can't show <= 0: pile on the edge
  d <- d[order(d$dens), ]                                # dense points drawn on top

  ybreaks <- seq(0, ymax, length.out = 4)
  ylabels <- if (is.null(ylab_max)) format(round(ybreaks), scientific = FALSE, trim = TRUE)
             else sprintf("%.1f", ybreaks * ylab_max / ymax)
  ylabels[1] <- "0"

  pt_layer <- if (rasterise_points)
    ggrastr::geom_point_rast(size = pt_size, shape = 16, raster.dpi = 300)
  else geom_point(size = pt_size, shape = 16)

  ggplot(d, aes(sybr, fsc, colour = dens)) +
    pt_layer +
    scale_colour_gradientn(colours = pal_ref, transform = "sqrt", guide = "none",
                           limits = c(1, max(cmax, 2)), oob = scales::squish) +
    scale_x_continuous(transform = "log10", limits = x_range, expand = c(0, 0),
                       breaks = 10^(0:ceiling(log10(x_range[2]))),
                       labels = scales::label_math(10^.x, format = log10),
                       name = "SYBR Green I") +
    scale_y_continuous(limits = c(0, ymax), oob = scales::squish, expand = c(0, 0),
                       breaks = ybreaks, labels = ylabels, name = "FSC") +
    annotate("rect", xmin = gate, xmax = x_range[2], ymin = 0, ymax = ymax,
             fill = NA, colour = "#2E7D32", linewidth = 0.6) +
    annotate("text", x = x_range[2], y = ymax, hjust = 1.1, vjust = 1.3, size = 3.2,
             lineheight = 0.95, colour = "grey15",
             label = sprintf("SYBR Green I\n%.2f%%", pct)) +
    labs(title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(size = 8, hjust = 0.5))
}

# Peak apexes: local maxima of a smoothed density on the log10 scale
find_peaks <- function(x, floor = peak_floor, adjust = peak_adjust, min_rel = peak_min_rel) {
  z <- log10(x[x > floor]); d <- density(z, adjust = adjust)
  i <- which(diff(sign(diff(d$y))) == -2) + 1
  out <- data.frame(position = 10^d$x[i], height = d$y[i])
  out[out$height >= min_rel * max(out$height), ]
}

# One peak, fitted in a window around its apex by maximum likelihood (no binning):
#   inside [lo, hi] (log10) the density is  p * truncated Normal(mu, s) + (1 - p) * baseline,
#   baseline flat, or sloped (slope = TRUE) so a neighbouring plateau isn't absorbed.
fit_peak <- function(x, center, w = 0.12, slope = FALSE) {
  z0 <- log10(center); lo <- z0 - w; hi <- z0 + w
  z  <- log10(x[x > 0]); z <- z[z > lo & z < hi]
  if (length(z) < 100) return(NULL)
  zc <- (lo + hi) / 2; half <- (hi - lo) / 2

  nll <- function(par) {
    mu <- par[1]; s <- exp(par[2]); p <- plogis(par[3])
    b  <- if (slope) tanh(par[4]) else 0                  # slope kept in (-1, 1): density stays > 0
    zt <- pnorm((hi - mu) / s) - pnorm((lo - mu) / s)     # truncation normaliser
    v  <- -sum(log(p * dnorm(z, mu, s) / zt +
                   (1 - p) * (1 + b * (z - zc) / half) / (hi - lo)))
    if (is.finite(v)) v else 1e10
  }
  start <- c(z0, log(if (slope) 0.08 else 0.05), 0, if (slope) 0)
  o  <- optim(start, nll, method = "Nelder-Mead", control = list(maxit = 3000))
  mu <- o$par[1]; s <- exp(o$par[2]); p <- plogis(o$par[3])
  b  <- if (slope) tanh(o$par[4]) else 0

  list(center = center, w = w, slope = slope, mu = mu, s = s, p = p, b = b, lo = lo, hi = hi,
       mean = 10^mu,
       cv_pct = sqrt(exp((s * log(10))^2) - 1) * 100,     # CV on the linear scale
       share_pct = 100 * p, n_window = length(z),
       # invalid if the Gaussian swallowed the window, left it, or the optimiser didn't converge
       valid = o$convergence == 0 && p < 0.99 && mu > lo && mu < hi)
}

# Try a narrow flat-baseline window first; if that fit isn't valid (broad peak on a plateau),
# retry with a wider window and a sloped baseline.
fit_peak_auto <- function(x, center) {
  for (sp in list(list(w = 0.12, slope = FALSE), list(w = 0.22, slope = TRUE))) {
    f <- fit_peak(x, center, sp$w, sp$slope)
    if (!is.null(f) && f$valid) return(f)
  }
  NULL
}

ci95 <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 20) c(NA_real_, NA_real_) else quantile(v, c(0.025, 0.975), names = FALSE)
}

# rbind data frames whose columns differ (missing columns filled with NA)
rbind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { for (cc in setdiff(cols, names(d))) d[[cc]] <- NA; d[cols] }))
}

# QC figure: histogram with each fitted Gaussian + baseline over its window
plot_fit_qc <- function(x, fits, title = NULL, bins = 120) {
  z   <- log10(x[x > 0])
  rng <- c(min(vapply(fits, `[[`, numeric(1), "lo")) - 0.15,
           max(vapply(fits, `[[`, numeric(1), "hi")) + 0.15)
  d   <- data.frame(z = z[z > rng[1] & z < rng[2]])
  curves <- do.call(rbind, lapply(seq_along(fits), function(i) {
    f  <- fits[[i]]; xs <- seq(f$lo, f$hi, length.out = 200)
    zt <- pnorm((f$hi - f$mu) / f$s) - pnorm((f$lo - f$mu) / f$s)
    zc <- (f$lo + f$hi) / 2; half <- (f$hi - f$lo) / 2
    sc <- f$n_window / nrow(d)                       # put curves on the histogram's density scale
    base <- (1 - f$p) * (1 + f$b * (xs - zc) / half) / (f$hi - f$lo)
    data.frame(x = xs, peak = factor(i),
               total = sc * (f$p * dnorm(xs, f$mu, f$s) / zt + base), baseline = sc * base)
  }))
  br <- seq(ceiling(rng[1]), floor(rng[2]))
  ggplot() +
    geom_histogram(data = d, aes(z, after_stat(density)), bins = bins, fill = "grey85", colour = NA) +
    geom_line(data = curves, aes(x, baseline, group = peak), colour = "grey40",
              linetype = "dashed", linewidth = 0.4) +
    geom_line(data = curves, aes(x, total, group = peak), colour = red, linewidth = 0.7) +
    scale_x_continuous(breaks = br, labels = parse(text = paste0("10^", br)), name = "SYBR Green I") +
    labs(y = "Density", title = title) +
    theme_classic(base_size = 11) + theme(plot.title = element_text(size = 8, hjust = 0.5))
}

## ---- 3. Read the files ------------------------------------------------------
message("== 1. Reading files ==")
files <- list.files(fcs_dir, pattern = file_pattern, full.names = TRUE,
                    recursive = TRUE, ignore.case = TRUE)
files <- files[order(tolower(basename(files)))]
if (length(files) == 0) stop("No files matching '", file_pattern, "' under '", fcs_dir, "'")
message(length(files), " file(s):\n  ", paste(basename(files), collapse = "\n  "))
ev <- lapply(files, read_events)
names(ev) <- basename(files)
is_control <- grepl(unstained_pattern, names(ev), ignore.case = TRUE)

## ---- 4. FACS panels ---------------------------------------------------------
message("== 2. FACS panels ==")
ymax <- max(vapply(ev, function(e) quantile(e$d$fsc, 0.995, names = FALSE), numeric(1)))
cmax <- quantile(unlist(lapply(ev, function(e) point_density(e$d, x_range, ymax, what = "cell"))),
                 0.99, names = FALSE)                  # shared colour scale across panels
panels <- lapply(seq_along(ev), function(i)
  plot_facs_panel(ev[[i]], title = tools::file_path_sans_ext(names(ev)[i]), ymax = ymax, cmax = cmax))
ncol_f <- min(3, length(panels))
facs_fig <- wrap_plots(panels, ncol = ncol_f) + plot_annotation(tag_levels = "A")
save_fig(file.path(out_dir, "facs_density.pdf"), facs_fig,
         width = 4 * ncol_f, height = 4 * ceiling(length(panels) / ncol_f))

facs_summary <- do.call(rbind, lapply(seq_along(ev), function(i) data.frame(
  file = names(ev)[i], control = is_control[i], n_events = nrow(ev[[i]]$d),
  n_above_gate = sum(ev[[i]]$d$sybr > gate), pct_above_gate = 100 * mean(ev[[i]]$d$sybr > gate))))
write.csv(facs_summary, file.path(out_dir, "facs_summary.csv"), row.names = FALSE)
cat("\nSYBR+ events (gate at ", gate, ")\n", sep = ""); print(facs_summary, digits = 3, row.names = FALSE)

## ---- 5. Peak analysis, one stained file at a time ----------------------------
analyse_peaks <- function(e) {
  fname <- basename(e$file)
  nm    <- gsub("[^A-Za-z0-9._-]+", "_", tools::file_path_sans_ext(fname))
  x     <- e$d$sybr; x <- x[x > gate]                     # events used for fitting and resampling
  message("\n== 3. Peaks: ", fname, " ==")

  centers <- if (!is.null(centers_override[[fname]])) centers_override[[fname]] else find_peaks(x)$position
  message("apexes: ", paste(round(centers, 1), collapse = ", "))
  fits <- Filter(Negate(is.null), lapply(centers, function(cn) fit_peak_auto(x, cn)))
  if (length(fits) == 0) { message("no peak could be fitted - skipping"); return(NULL) }
  fits <- fits[order(vapply(fits, `[[`, numeric(1), "mean"))]
  K <- length(fits); means <- vapply(fits, `[[`, numeric(1), "mean")
  if (K < length(centers)) message((length(centers) - K), " apex(es) dropped: no valid fit")

  ## bootstrap: resample events, refit every peak with the SAME window/baseline, repeat
  message("== 4. Bootstrap (", n_boot, " resamples) ==")
  set.seed(seed)
  boot <- matrix(NA_real_, n_boot, K)
  for (i in seq_len(n_boot)) {
    xb <- sample(x, replace = TRUE)
    for (k in seq_len(K)) {
      f <- tryCatch(fit_peak(xb, fits[[k]]$center, fits[[k]]$w, fits[[k]]$slope),
                    error = function(err) NULL)
      if (!is.null(f) && f$valid) boot[i, k] <- f$mean
    }
    if (i %% 100 == 0) message("  ", i, "/", n_boot)
  }

  tab_pos <- data.frame(
    file = fname, peak = seq_len(K),
    apex = vapply(fits, `[[`, numeric(1), "center"), mean = means,
    ci_low  = vapply(seq_len(K), function(k) ci95(boot[, k])[1], numeric(1)),
    ci_high = vapply(seq_len(K), function(k) ci95(boot[, k])[2], numeric(1)),
    cv_pct = vapply(fits, `[[`, numeric(1), "cv_pct"),
    share_pct = vapply(fits, `[[`, numeric(1), "share_pct"),
    window_half_width_log10 = vapply(fits, `[[`, numeric(1), "w"),
    baseline = ifelse(vapply(fits, `[[`, logical(1), "slope"), "sloped", "flat"),
    valid_fits = sprintf("%d/%d", colSums(!is.na(boot)), n_boot))

  tab_fold <- NULL
  if (K > 1) tab_fold <- do.call(rbind, lapply(seq_len(K - 1), function(i) {
    rb <- boot[, i + 1] / boot[, i]                       # resamples where BOTH peaks fitted validly
    r  <- means[i + 1] / means[i]; ci <- ci95(rb)
    data.frame(file = fname, from = i, to = i + 1, fold_change = r,
               log10_distance = log10(r), log2_distance = log2(r),   # exact doubling: 0.301 / 1.000
               ci_low = ci[1], ci_high = ci[2], valid_fits = sprintf("%d/%d", sum(is.finite(rb)), n_boot),
               ci_includes_2x = if (is.na(ci[1])) NA else ci[1] <= 2 & ci[2] >= 2)
  }))

  ## sensitivity to the window choice (the bootstrap does not cover this)
  tab_sens <- do.call(rbind, lapply(c(0.8, 1, 1.25), function(a) {
    m <- vapply(fits, function(f) {
      g <- fit_peak(x, f$center, f$w * a, f$slope)
      if (is.null(g) || !g$valid) NA_real_ else g$mean }, numeric(1))
    out <- data.frame(file = fname, window_scale = a, t(setNames(m, paste0("peak", seq_len(K)))))
    if (K > 1) for (i in seq_len(K - 1)) out[[paste0("fold_", i, "_", i + 1)]] <- m[i + 1] / m[i]
    out
  }))

  cat("\nPeak positions (95% bootstrap interval)\n"); print(tab_pos, digits = 3, row.names = FALSE)
  if (!is.null(tab_fold)) { cat("\nFold change between neighbouring peaks\n"); print(tab_fold, digits = 3, row.names = FALSE) }
  cat("\nSensitivity to window width\n"); print(tab_sens, digits = 3, row.names = FALSE)

  ## figures
  save_fig(file.path(out_dir, paste0(nm, "_peak_fit_qc.pdf")),
           plot_fit_qc(x, fits, title = fname), width = 6, height = 3.5)

  # density (top) above the FACS panel (bottom), same log10 x axis, fitted peak means through both
  xd <- pmin(x, x_range[2])
  pd <- ggplot(data.frame(x = xd), aes(x)) +
    geom_density(aes(y = after_stat(scaled)), fill = "grey85", colour = "grey20",
                 adjust = 0.6, linewidth = 0.5) +
    geom_vline(xintercept = means, linetype = "dotted", colour = red) +
    scale_x_log10(limits = x_range, expand = c(0, 0), breaks = 10^(0:ceiling(log10(x_range[2]))),
                  labels = scales::label_math(10^.x, format = log10)) +
    scale_y_continuous(breaks = c(0, 0.5, 1), name = "Rel. density") +
    labs(title = tools::file_path_sans_ext(fname)) +
    theme_classic(base_size = 11) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(size = 8, hjust = 0.5))
  cap <- NULL
  if (K > 1) {
    dd <- data.frame(from = head(means, -1), to = tail(means, -1)); dd$fold <- dd$to / dd$from
    pd <- pd +
      geom_segment(data = dd, aes(x = from, xend = to, y = 1.10, yend = 1.10), inherit.aes = FALSE,
                   colour = "grey15", arrow = arrow(ends = "both", length = unit(0.05, "in"))) +
      geom_text(data = dd, aes(x = sqrt(from * to), y = 1.15, label = sprintf("%.2f\u00d7", fold)),
                inherit.aes = FALSE, vjust = 0, size = 2.8)
    cap <- paste0("Fold change, 95% bootstrap interval (", n_boot, " resamples):\n",
                  paste(sprintf("peak %d to peak %d: %.2f\u00d7 [%.2f\u2013%.2f]",
                                tab_fold$from, tab_fold$to, tab_fold$fold_change,
                                tab_fold$ci_low, tab_fold$ci_high), collapse = "\n"))
  }
  pd <- pd + coord_cartesian(ylim = c(0, 1.4), expand = FALSE)
  pf <- plot_facs_panel(e, title = NULL, ymax = ymax, cmax = cmax) +
    geom_vline(xintercept = means, linetype = "dotted", colour = red)
  stack <- (pd / pf) + plot_layout(heights = c(1, 2.6)) +
    plot_annotation(caption = cap, theme = theme(plot.caption = element_text(size = 8, hjust = 0)))
  save_fig(file.path(out_dir, paste0(nm, "_density_facs_stack.pdf")), stack, width = 5, height = 7)

  list(pos = tab_pos, fold = tab_fold, sens = tab_sens)
}

res <- lapply(ev[!is_control], analyse_peaks)
res <- Filter(Negate(is.null), res)

## ---- 6. Tables and session info ----------------------------------------------
if (length(res)) {
  write.csv(do.call(rbind, lapply(res, `[[`, "pos")),  file.path(out_dir, "peak_positions.csv"),   row.names = FALSE)
  folds <- Filter(Negate(is.null), lapply(res, `[[`, "fold"))
  if (length(folds)) write.csv(do.call(rbind, folds),  file.path(out_dir, "peak_fold_changes.csv"), row.names = FALSE)
  write.csv(rbind_fill(lapply(res, `[[`, "sens")),   file.path(out_dir, "window_sensitivity.csv"), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
message("\nDone. Outputs in '", out_dir, "/':\n  ", paste(list.files(out_dir), collapse = "\n  "))
