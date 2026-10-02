# Radar Target Generation and Detection
# R port of src/radar-target-generation-and-detection.m
#
# This file is a port of the MATLAB template.
#
# Run:   Rscript src/radar-target-generation-and-detection.R
# Needs: base R only (stats + graphics).

## Setup
suppressWarnings(dir.create("plots", showWarnings = FALSE))
save_png <- !interactive()
# Plotting helpers
open_device <- function(name, width = 1000, height = 800) {
  if (save_png) {
    png(file.path("plots", paste0(name, ".png")),
        width = width, height = height)
  } else {
    grDevices::dev.new()
  }
}
close_device <- function() {
  if (save_png) invisible(grDevices::dev.off())
}
# Facet colors for persp surface plots: map each facet's height (the
# average of its four corner cells) linearly onto a palette, like a
# heatmap. persp expects one color per facet, i.e. (nrow-1)*(ncol-1).
# The color scale is derived from the facet heights themselves — not from
# the cell values — because a detection cluster only one cell wide never
# produces a facet of full cell height (its facets average to 0.5 and
# would otherwise stay below the top color).
facet_cols <- function(z, palette) {
  nr <- nrow(z); nc <- ncol(z)
  h <- (z[-nr, -nc] + z[-nr, -1] + z[-1, -nc] + z[-1, -1]) / 4
  if (max(h) == min(h)) return(rep(palette[1], (nr - 1) * (nc - 1)))
  breaks <- seq(min(h), max(h), length.out = length(palette) + 1)
  palette[cut(h, breaks = breaks, include.lowest = TRUE)]
}
# Vertical color-scale legend for a facet heatmap: a thin strip of the
# palette at the right edge of the figure, labelled with the mapped values
# (bottom = min(z), top = max(z)). Call after the persp plot, with the same
# z and palette.
add_color_legend <- function(z, palette, label = "", at = NULL) {
  lo <- min(z); hi <- max(z)
  n <- length(palette)
  if (is.null(at)) at <- pretty(c(lo, hi))
  at <- at[at >= lo & at <= hi]
  op <- par(fig = c(0.90, 0.97, 0.12, 0.88), new = TRUE,
            mar = c(0.3, 0.1, 0.3, 2.2))
  on.exit(par(op), add = TRUE)
  graphics::image(1, seq_len(n), matrix(seq_len(n), nrow = 1),
                 col = palette, axes = FALSE, xlab = "", ylab = "")
  graphics::box(col = "grey50")
  graphics::axis(4, at = 1 + (at - lo) / (hi - lo) * (n - 1), labels = at,
                 cex.axis = 0.7, tick = FALSE)
  graphics::mtext(label, side = 4, line = 1.6, cex = 0.7)
}

## Radar Specifications
############################
# Frequency of operation = 77GHz
# Max Range = 200m
# Range Resolution = 1 m
# Max Velocity = 100 m/s
############################
c          <- 3e8     # speed of light (m/s)
fc         <- 77e9    # operating carrier frequency (Hz)
max_range  <- 200     # maximum range (m)
range_res  <- 1       # range resolution (m)
max_vel    <- 100     # maximum velocity (m/s)

## User Defined Range and Velocity of target
# Initial position and constant radial velocity of the target
# (negative = approaching).
target_range0 <- 75    # initial range (m)
target_vel    <- -25   # constant radial velocity (m/s)

## FMCW Waveform Generation
# Design the FMCW waveform from the requirements above: Bandwidth (B),
# Chirp Time (Tchirp) and Slope (slope) of the FMCW chirp.
B      <- c / (2 * range_res)             # bandwidth (Hz): 150 MHz
Tchirp <- 5.5 * (2 * max_range / c)       # chirp time (s)
slope  <- B / Tchirp                      # chirp slope (Hz/s)

# Report the designed waveform parameters.
cat("FMCW waveform design:\n")
cat(sprintf("  Bandwidth  B      = %.1f MHz\n", B / 1e6))
cat(sprintf("  Chirp time Tchirp = %.3f us\n", Tchirp * 1e6))
cat(sprintf("  Slope      slope  = %.4e Hz/s\n", slope))

# Operating carrier frequency of Radar
fc <- 77e9             # carrier freq

# The number of chirps in one sequence. Its ideal to have 2^ value for the ease
# of running the FFT for Doppler Estimation.
Nd <- 128              # #of doppler cells OR #of sent periods  # number of chirps

# The number of samples on each chirp.
Nr <- 1024             # for length of time OR # of range cells

# Timestamp for running the displacement scenario for every sample on each
# chirp.
t <- seq(0, Nd * Tchirp, length.out = Nr * Nd)   # total time for samples

# Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx <- numeric(length(t))   # transmitted signal
Rx <- numeric(length(t))   # received signal
Mix <- numeric(length(t))  # beat signal

# Similar vectors for range_covered and time delay.
r_t <- numeric(length(t))
td  <- numeric(length(t))

## FMCW Waveform Visualization
# Plots the previously designed FMCW waveform (B, Tchirp, slope from above)
# on the previously defined timestamp vector t: the transmitted chirp ramp,
# the delayed receive ramp for the target's initial position, and the
# resulting (nearly constant) beat frequency that the range FFT will later
# measure. The actual Tx/Rx/Mix signals are generated in the simulation
# loop below.

td0 <- 2 * target_range0 / c          # round-trip delay at initial position (s)
fb  <- slope * td0                    # beat frequency from range (Hz)

open_device("00_fmcw_chirp", width = 1000, height = 1200)
par(mfrow = c(3, 1), mar = c(5, 5, 3, 2) + 0.1)

# (1) Transmitted chirp: instantaneous frequency fc + slope*t over two
# chirps (first 2*Nr samples of the defined timestamp vector t).
t2 <- t[1:(2 * Nr)]
plot(t2 * 1e6, (fc + slope * (t2 %% Tchirp)) / 1e9, type = "l", col = "blue",
     lwd = 2, main = "FMCW Waveform: Transmitted Chirp",
     xlab = "Time [us]", ylab = "Frequency [GHz]")
abline(h = c(fc, fc + B) / 1e9, col = "gray", lty = 3)
legend("topleft",
       sprintf("slope = %.3e Hz/s (B = %.0f MHz, Tchirp = %.2f us)",
               slope, B / 1e6, Tchirp * 1e6), bty = "n")

# (2) Tx ramp vs delayed Rx ramp within one chirp (first Nr samples of t):
# the vertical offset between the ramps is the beat frequency
# fb = slope * td. The Doppler shift of the moving target is omitted here —
# it does not change the range beat, it shows up as a phase shift between
# consecutive chirps.
t1 <- t[1:Nr]
plot(t1 * 1e6, (fc + slope * t1) / 1e9, type = "l", col = "blue", lwd = 2,
     main = "FMCW Waveform: Tx vs Rx (delayed) Chirp",
     xlab = "Time [us]", ylab = "Frequency [GHz]")
lines(t1 * 1e6, (fc + slope * (t1 - td0)) / 1e9, col = "red", lwd = 2)
abline(v = td0 * 1e6, col = "gray", lty = 3)
x_mark <- 0.85 * max(t1)
arrows(x_mark * 1e6, (fc + slope * x_mark) / 1e9,
       x_mark * 1e6, (fc + slope * (x_mark - td0)) / 1e9,
       code = 3, angle = 90, length = 0.05)
text(x_mark * 1e6, (fc + slope * (x_mark - td0 / 2)) / 1e9,
     sprintf("fb = %.2f MHz -> R = %.1f m", fb / 1e6, fb * c / (2 * slope)),
     pos = 2, cex = 0.9, offset = 1)
legend("topleft", c("Tx", "Rx (delayed by td = 2R/c)"),
       col = c("blue", "red"), lty = 1, lwd = 2, bty = "n")

# (3) Beat frequency over one chirp: essentially constant, it encodes range.
# The approaching target lets fb fall by only ~25 Hz over the chirp
# (vs ~10.2 MHz), so the y axis is fixed to the full scale: the curve is
# flat for all practical purposes, and the range FFT sees one sharp bin.
fb_t <- slope * 2 * (target_range0 + target_vel * t1) / c
plot(t1 * 1e6, fb_t / 1e6, type = "l", col = "darkgreen", lwd = 2,
     main = "FMCW Waveform: Beat Frequency (slope * 2R(t)/c)",
     xlab = "Time [us]", ylab = "Beat frequency [MHz]",
     ylim = c(0, 1.05 * max(fb_t) / 1e6))
legend("bottomleft",
       sprintf("fb = %.2f MHz -> R = fb*c/(2*slope) = %.1f m\n(falls by only %.0f Hz over the chirp as the target approaches)",
               fb / 1e6, fb * c / (2 * slope),
               slope * 2 * abs(target_vel) * Tchirp / c), bty = "n")

close_device()

## Signal generation and Moving Target simulation
# Running the radar scenario over the time.

# Remark: The below code is vectorized to better show the simplicity of the signal formation

# Range of the target for constant velocity, at every time stamp.
r_t <- target_range0 + target_vel * t
# Time delay of the received signal (out and back).
td <- 2 * r_t / c
# Theoretical Doppler shift of the moving target (course formula
# fd = 2*fc*v/c) — kept for the comparison with the simulation.
fd <- 2 * target_vel * fc / c

# Signal path per wave equations of the tutorial:
# send the FMCW chirp,
# receive the time-delayed reflection (t replaced by t - td),
# and mix both.
# The element-wise product Tx * Rx implements the
# mixing of the two signals ("frequency subtraction").
#   Tx = cos(2*pi*(fc*t + 0.5*slope*t^2))
#   Rx = cos(2*pi*(fc*(t-td) + 0.5*slope*(t-td)^2))
Tx <- cos(2 * pi * (fc * t + 0.5 * slope * t^2))
Rx <- cos(2 * pi * (fc * (t - td) + 0.5 * slope * (t - td)^2))

# Beat signal: element-wise product of transmitted and received signal.
Mix <- Tx * Rx

# Visualize the generated beat signal in the time domain:
# the mixing product oscillates dominantly at fb + fd — the range beat frequency
# plus the small Doppler shift (10.21 MHz for the target's initial
# position) — superimposed with the fast sum-frequency term of the analog
# product. Its period is the range information that the range FFT will
# extract as a peak in task 3.

open_device("00_beat_signal", width = 1000, height = 1000)
par(mfrow = c(2, 1), mar = c(5, 5, 3, 2) + 0.1)

# (1) Complete first chirp: the beat signal oscillates at fb + fd.
plot(t1 * 1e6, Mix[1:Nr], type = "l", col = "darkblue",
     main = "Beat Signal (Mix = Tx * Rx), First Chirp",
     xlab = "Time [us]", ylab = "Amplitude")
abline(h = 0, col = "gray", lty = 3)
legend("topright",
       sprintf("beat tone at fb + fd = %.2f MHz (R = %.1f m)",
               (fb + fd) / 1e6, fb * c / (2 * slope)), bty = "n")

# (2) Zoom on the first microsecond: the individual beat cycles.
n_zoom <- which(t1 >= 1e-6)[1]
plot(t1[1:n_zoom] * 1e6, Mix[1:n_zoom], type = "l", col = "darkblue",
     main = "Beat Signal, First 1 us (zoom)",
     xlab = "Time [us]", ylab = "Amplitude")
abline(h = 0, col = "gray", lty = 3)
legend("bottomright",
       sprintf("period ~ %.0f ns -> fb + fd ~ %.2f MHz -> R = %.1f m",
               1 / (fb + fd) * 1e9, (fb + fd) / 1e6, fb * c / (2 * slope)),
       bty = "n")

close_device()


## RANGE MEASUREMENT

# Reshape the beat signal into an Nr x Nd array (column-major, matching
# MATLAB's reshape: column j holds chirp j).
Mix2D <- matrix(Mix, nrow = Nr, ncol = Nd)

# FFT of the first chirp along the range dimension, normalized by Nr.
sig_fft <- fft(Mix2D[, 1]) / Nr
# Absolute value of the FFT output.
sig_fft <- Mod(sig_fft)
# Output of FFT is double sided; keep only one side of the spectrum.
sig_fft <- sig_fft[1:(Nr / 2)]

range_fft_axis <- (0:(Nr / 2 - 1)) * (c / (2 * B))   # bin k -> k * range resolution (m)

# Report the detected range peak.
k_peak <- which.max(sig_fft) - 1
cat(sprintf("Range FFT peak: bin %d -> %.1f m (target_range0 = %.0f m)\n",
            k_peak, range_fft_axis[k_peak + 1], target_range0))

open_device("01_range_fft")
plot(range_fft_axis, sig_fft / max(sig_fft), type = "l",
     xlim = c(0, 200), ylim = c(0, 1),
     xlab = "Range (m)", ylab = "Normalized amplitude",
     main = "Range from First FFT")
grid(col = "grey80")
close_device()

## RANGE DOPPLER RESPONSE
# 2D FFT on the beat signal generates the Range-Doppler Map (RDM).

# 2D FFT: FFT along the range dimension (columns), then along Doppler (rows).
sig_fft2 <- t(apply(apply(Mix2D, 2, fft), 1, fft))
# One side of the spectrum in the range dimension, then fftshift (both dims).
sig_fft2 <- sig_fft2[1:(Nr / 2), ]
nr_half  <- Nr / 2
nc_half  <- Nd %/% 2
sig_fft2 <- sig_fft2[c((nr_half %/% 2 + 1):nr_half, 1:(nr_half %/% 2)),
                     c((nc_half + 1):Nd, 1:nc_half)]  # == MATLAB fftshift
RDM <- 10 * log10(Mod(sig_fft2))
# clamp -Inf (zero-power bins) so plotting works
RDM[!is.finite(RDM)] <- min(RDM[is.finite(RDM)])

# Axis conversion from bins to range and velocity based on their max values.
doppler_axis <- seq(-100, 100, length.out = Nd)
range_axis   <- seq(-200, 200, length.out = Nr / 2) * ((Nr / 2) / 400)

# Convert an RDM cell (row/column after fftshift) back to physical range and
# velocity:
# undo the fftshift to recover the original FFT bins,
# then range = bin * c/(2B) (1 m per range bin) and fd = bin/(Nd*Tchirp),
# v = fd*c/(2*fc).
cell_to_range_velocity <- function(cell) {
  r_bin <- if (cell[1] <= Nr / 4) Nr / 4 + cell[1] else cell[1] - Nr / 4
  d_bin <- if (cell[2] <= nc_half) nc_half + cell[2] else cell[2] - nc_half
  fd <- if (d_bin - 1 <= Nd / 2) (d_bin - 1) / (Nd * Tchirp) else
    (d_bin - 1 - Nd) / (Nd * Tchirp)
  list(range = (r_bin - 1) * c / (2 * B), velocity = fd * c / (2 * fc))
}

# Report the strongest RDM cell: the target peak in physical units.
tgt <- cell_to_range_velocity(arrayInd(which.max(RDM), dim(RDM)))
cat(sprintf("RDM peak (target): range = %.1f m, velocity = %.1f m/s\n",
            tgt$range, tgt$velocity))
# Sidenote: the theoretical formula of the course, fd = 2*fc*v/c, predicts
# the Doppler velocity directly. The measured Doppler deviates from it
# because the theory holds the range constant while the simulation
# updates the range at every time stamp (see the report, section 4.1).
cat(sprintf("Theoretical velocity from fd = 2*fc*v/c: %.1f m/s\n",
            fd * c / (2 * fc)))

open_device("02_range_doppler_map")
graphics::filled.contour(doppler_axis, range_axis, t(RDM),
                         color.palette = grDevices::hcl.colors,
                         xlab = "Velocity (m/s)", ylab = "Range (m)",
                         main = "Range Doppler Map (2D FFT)")
close_device()

open_device("03_range_doppler_surface")
rdm_palette <- grDevices::hcl.colors(64, "YlOrRd")
par(mar = c(5, 4, 4, 7))   # reserve the right margin for the color scale
graphics::persp(doppler_axis, range_axis, t(RDM),
                theta = 25, phi = 15, shade = NA,
                col = facet_cols(t(RDM), rdm_palette),
                border = NA, ticktype = "detailed",
                xlab = "Velocity (m/s)", ylab = "Range (m)", zlab = "dB",
                main = "Range Doppler Response")
add_color_legend(t(RDM), rdm_palette, label = "dB")
close_device()

## CFAR implementation
# Slide a window through the complete Range-Doppler Map.

# Number of Training Cells in both dimensions.
Tr <- 10   # range training cells (per side)
Td <- 4    # Doppler training cells (per side)
# Number of Guard Cells around the Cell Under Test (CUT).
Gr <- 3    # range guard cells (per side)
Gd <- 20    # Doppler guard cells (per side)
# Threshold offset by SNR value in dB (14 dB lets the detection cluster
# cover the full Doppler smear of the target — including its leading
# edge, which carries the frame-start Doppler — without picking up
# leakage cells; see the report for the measured offset sweep).
offset <- 14

# Linear power version of the map: db2pow(x) = 10^(x/10).
RDM_pow <- 10^(RDM / 10)

# CFAR output, same size as the RDM (untested edge cells stay 0).
CFAR <- matrix(0, nrow = Nr / 2, ncol = Nd)

# The training cells form the ring between the training window and the
# guard block around the CUT: (2*(Tr+Gr)+1) x (2*(Td+Gd)+1) cells minus the
# (2*Gr+1) x (2*Gd+1) guard block (which contains the CUT itself).
n_train <- (2 * (Tr + Gr) + 1) * (2 * (Td + Gd) + 1) -
           (2 * Gr + 1) * (2 * Gd + 1)

# Slide the CUT across the RDM, leaving margins of Tr+Gr cells in the range
# dimension and Td+Gd cells in the Doppler dimension, so that the full
# training window always fits inside the map.
for (i in (Tr + Gr + 1):(Nr / 2 - Tr - Gr)) {
  for (j in (Td + Gd + 1):(Nd - Td - Gd)) {
    # Average noise power over the training ring: sum over the full
    # training window, subtract the guard block, average by the number of
    # training cells.
    noise_level <- (sum(RDM_pow[(i - Tr - Gr):(i + Tr + Gr),
                                (j - Td - Gd):(j + Td + Gd)]) -
                    sum(RDM_pow[(i - Gr):(i + Gr), (j - Gd):(j + Gd)])) /
                   n_train
    # Threshold: averaged noise power back to dB (pow2db) plus the offset.
    threshold_db <- 10 * log10(noise_level) + offset
    # CUT above threshold -> 1, else 0.
    CFAR[i, j] <- as.numeric(RDM[i, j] > threshold_db)
  }
}

## CFAR output display (surf / image like the Range Doppler Response)

open_device("04_cfar_image")
graphics::image(doppler_axis, range_axis, t(CFAR), col = c("white", "red"),
                xlab = "Velocity (m/s)", ylab = "Range (m)",
                main = sprintf("CFAR Detection Output (%d cells)", sum(CFAR)))
graphics::legend("topright", legend = c("no detection", "detection"),
                 fill = c("white", "red"), bty = "n")
close_device()

open_device("05_cfar_surf")
cfar_palette <- c("grey90", "red")
graphics::persp(doppler_axis, range_axis, t(CFAR),
                theta = 25, phi = 15, shade = NA,
                col = facet_cols(t(CFAR), cfar_palette),
                border = "grey40", ticktype = "detailed",
                xlab = "Velocity (m/s)", ylab = "Range (m)", zlab = "detected",
                main = "CFAR Output")
add_color_legend(t(CFAR), cfar_palette, label = "detected", at = c(0, 1))
close_device()

## Sanity summary on the console.
cat(sprintf("CFAR detections: %d of %d cells\n", sum(CFAR), length(CFAR)))

# Strongest detection: the detected cell with the highest RDM level,
# reported in physical units.
det_cells <- which(CFAR == 1, arr.ind = TRUE)
if (nrow(det_cells) > 0) {
  best <- det_cells[which.max(RDM[det_cells]), ]
  best_phys <- cell_to_range_velocity(best)
  cat(sprintf("Strongest detection: range = %.1f m, velocity = %.1f m/s\n",
              best_phys$range, best_phys$velocity))

  # Leading edge of the Doppler smear: within the frame the moving
  # target's beat frequency drifts steadily in the direction of its
  # relative velocity (toward more negative Doppler for the approaching
  # target here), so the per-chirp phase increment at the frame start is
  # the pure Doppler shift fd = 2*fc*v/c. The edge of the detection
  # cluster on the undrifted side therefore recovers the target's
  # undrifted velocity — the least-negative detected Doppler bin for the
  # approaching target.
  det_vel <- apply(det_cells, 1, function(cell) {
    cell_to_range_velocity(cell)$velocity
  })
  cat(sprintf("Doppler leading edge: %.1f m/s (frame-start Doppler)\n",
              max(det_vel)))

  # Values of the selected (detected) target bins: range, velocity and
  # level of every detection, sorted by velocity. The listing shows that
  # the detected bins cover the target: all bins sit at the target's
  # range bin (75 m) and span the smeared Doppler up to the leading
  # edge, which carries the target's velocity.
  cat("Detected target bins:\n")
  for (k in order(det_vel)) {
    p <- cell_to_range_velocity(det_cells[k, ])
    cat(sprintf("  range = %5.1f m, velocity = %6.1f m/s, level = %4.1f dB\n",
                p$range, p$velocity, RDM[det_cells[k, 1], det_cells[k, 2]]))
  }
}
