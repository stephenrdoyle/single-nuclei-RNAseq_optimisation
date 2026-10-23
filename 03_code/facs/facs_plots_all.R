## =============================================================================
## facs_plots_all.R
##
## FACS density plots (FSC vs SYBR Green I) for every .fcs file in a folder.
##
## Writes to facs_plots/:
##   facs_density_all.pdf   all files in one multi-page PDF (3 columns x 4 rows per page)
##   individual/<file>.pdf  one PDF per file (set individual_pdfs <- FALSE to skip)
##   facs_summary.csv       events and % above the SYBR+ gate for every file
##
## Usage (from the project folder):   source("facs_plots_all.R")
##                                or: Rscript facs_plots_all.R
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
fcs_dir           <- "raw_fcs"        # searched recursively
file_pattern      <- "\\.fcs$"        # regex on file names; e.g. "tmuris.*\\.fcs$" for a subset
unstained_pattern <- "unstained"      # only used to label controls in the summary table
sybr_ch           <- "FL53-A"         # SYBR Green I, 507/19 detector
fsc_ch            <- "FSC04-A"        # 488 nm forward scatter
out_dir           <- "facs_plots"

gate      <- 10            # SYBR+ gate (green box); the percentage above it is printed on each panel
x_range   <- c(1, 1e5)     # SYBR axis, log10: 10^0 - 10^5
ylab_max  <- 512           # relabel FSC ticks to a 0-512 display scale (labels only); NULL = raw values
shared_axes <- TRUE        # TRUE: same FSC range and colour scale on every panel (comparable);
                           # FALSE: each panel scaled to its own data
individual_pdfs  <- TRUE
rasterise_points <- FALSE  # TRUE: dots as an image inside the PDF (needs ggrastr), smaller files
ncol_page <- 3; nrow_page <- 4          # panels per page of the combined PDF

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
    stop("channels ", sybr_ch, "/", fsc_ch, " not found. Available: ", paste(cn, collapse = ", "))
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

## ---- 3. Read the files ------------------------------------------------------
message("== Reading files ==")
files <- list.files(fcs_dir, pattern = file_pattern, full.names = TRUE,
                    recursive = TRUE, ignore.case = TRUE)
files <- files[order(tolower(basename(files)))]
if (length(files) == 0) stop("No files matching '", file_pattern, "' under '", fcs_dir, "'")
message(length(files), " file(s) found")

keys <- make.unique(basename(files))           # unique names even if two folders share a file name
ev <- list()
for (i in seq_along(files)) {
  r <- tryCatch(read_events(files[i]),
                error = function(err) { message("  SKIPPED ", keys[i], ": ", conditionMessage(err)); NULL })
  if (!is.null(r)) ev[[keys[i]]] <- r
}
if (length(ev) == 0) stop("No file could be read.")
message(length(ev), " file(s) read, ", length(files) - length(ev), " skipped")

## ---- 4. Scales shared across panels -------------------------------------------
if (shared_axes) {
  ymax_all <- max(vapply(ev, function(e) quantile(e$d$fsc, 0.995, names = FALSE), numeric(1)))
  cmax_all <- quantile(unlist(lapply(ev, function(e)
    point_density(e$d, x_range, ymax_all, what = "cell"))), 0.99, names = FALSE)
  message("Shared scales: FSC 0 to ", signif(ymax_all, 3), "; colour cap ", signif(cmax_all, 3),
          " events per cell")
}

## ---- 5. Panels --------------------------------------------------------------
message("== Plotting ==")
wrap_title <- function(s) paste(strwrap(s, 42), collapse = "\n")
panels <- list()
for (k in names(ev)) {
  e <- ev[[k]]
  ymax <- if (shared_axes) ymax_all else quantile(e$d$fsc, 0.995, names = FALSE)
  cmax <- if (shared_axes) cmax_all else
    quantile(point_density(e$d, x_range, ymax, what = "cell"), 0.99, names = FALSE)
  panels[[k]] <- plot_facs_panel(e, title = wrap_title(tools::file_path_sans_ext(k)),
                                 ymax = ymax, cmax = cmax)
}

# one multi-page PDF with all panels
per_page <- ncol_page * nrow_page
pages <- split(seq_along(panels), ceiling(seq_along(panels) / per_page))
all_pdf <- file.path(out_dir, "facs_density_all.pdf")
if (isTRUE(capabilities("cairo"))) {
  grDevices::cairo_pdf(all_pdf, width = 4 * ncol_page, height = 4 * nrow_page, onefile = TRUE)
} else {
  grDevices::pdf(all_pdf, width = 4 * ncol_page, height = 4 * nrow_page, onefile = TRUE)
}
for (idx in pages)
  print(wrap_plots(panels[idx], ncol = ncol_page, nrow = nrow_page) +
          plot_annotation(tag_levels = "A"))
invisible(dev.off())
message("Saved ", all_pdf, " (", length(pages), " page(s))")

# one PDF per file
if (individual_pdfs) {
  dir.create(file.path(out_dir, "individual"), showWarnings = FALSE)
  for (k in names(panels))
    save_fig(file.path(out_dir, "individual",
                       paste0(gsub("[^A-Za-z0-9._-]+", "_", tools::file_path_sans_ext(k)), ".pdf")),
             panels[[k]], width = 4, height = 4)
  message("Saved ", length(panels), " individual PDFs in ", file.path(out_dir, "individual"))
}

## ---- 6. Summary table -------------------------------------------------------
summ <- do.call(rbind, lapply(names(ev), function(k) data.frame(
  file = k, control = grepl(unstained_pattern, k, ignore.case = TRUE),
  n_events = nrow(ev[[k]]$d), n_above_gate = sum(ev[[k]]$d$sybr > gate),
  pct_above_gate = 100 * mean(ev[[k]]$d$sybr > gate))))
write.csv(summ, file.path(out_dir, "facs_summary.csv"), row.names = FALSE)
cat("\nSYBR+ events (gate at ", gate, ")\n", sep = ""); print(summ, digits = 3, row.names = FALSE)
message("\nDone. Outputs in '", out_dir, "/'")
