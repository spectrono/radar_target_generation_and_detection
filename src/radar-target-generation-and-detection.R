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

# *TODO*:
# define the target's initial position and velocity. Note : Velocity
# remains constant.
target_range0 <- 75    # initial range (m)
target_vel    <- -25   # constant radial velocity (m/s), negative = approaching

## FMCW Waveform Generation

# *TODO* :
# Design the FMCW waveform by giving the specs of each of its parameters.
# Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
# chirp using the requirements above.
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

# Transmitted and received signal for every time sample.
Tx <- cos(2 * pi * (fc * t + 0.5 * slope * t^2))
Rx <- cos(2 * pi * (fc * (t - td) + 0.5 * slope * (t - td)^2))

# Beat signal: element-wise product of Transmit and Receive signal.
Mix <- Tx * Rx

# Visualize the generated beat signal in the time domain: the product
# Tx * Rx oscillates dominantly at the beat frequency fb = 10.23 MHz
# (superimposed with the fast sum-frequency term of the mixing product).
# Its period is the range information that the range FFT will extract as a
# peak in task 3.

open_device("00_beat_signal", width = 1000, height = 1000)
par(mfrow = c(2, 1), mar = c(5, 5, 3, 2) + 0.1)

# (1) Complete first chirp: the beat signal oscillates at fb.
plot(t1 * 1e6, Mix[1:Nr], type = "l", col = "darkblue",
     main = "Beat Signal (Mix = Tx * Rx), First Chirp",
     xlab = "Time [us]", ylab = "Amplitude")
abline(h = 0, col = "gray", lty = 3)
legend("topright",
       sprintf("dominant oscillation at fb = %.2f MHz (R = %.1f m)",
               fb / 1e6, fb * c / (2 * slope)), bty = "n")

# (2) Zoom on the first microsecond: the individual beat cycles.
n_zoom <- which(t1 >= 1e-6)[1]
plot(t1[1:n_zoom] * 1e6, Mix[1:n_zoom], type = "l", col = "darkblue",
     main = "Beat Signal, First 1 us (zoom)",
     xlab = "Time [us]", ylab = "Amplitude")
abline(h = 0, col = "gray", lty = 3)
legend("bottomright",
       sprintf("period ~ %.0f ns -> fb ~ %.2f MHz -> R = %.1f m",
               1 / fb * 1e9, fb / 1e6, fb * c / (2 * slope)), bty = "n")

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
# The 2D FFT implementation is already provided here. This will run a 2DFFT
# on the mixed signal (beat signal) output and generate a range doppler
# map. You will implement CFAR on the generated RDM.

# Range Doppler Map Generation.

# The output of the 2D FFT is an image that has response in range and
# doppler FFT bins. So, it is important to convert the axis from bin sizes
# to range and doppler based on their Max values.

Mix2D <- matrix(Mix, nrow = Nr, ncol = Nd)

# 2D FFT using the FFT size for both dimensions:
# FFT along the range dimension (columns), then along Doppler (rows).
sig_fft2 <- t(apply(apply(Mix2D, 2, fft), 1, fft))

# Taking just one side of signal from Range dimension, then fftshift
# (circular shift by half the length, both dimensions == MATLAB fftshift).
nr_half <- Nr / 2
nc_half <- Nd %/% 2
sig_fft2 <- sig_fft2[1:(Nr / 2), ]
sig_fft2 <- sig_fft2[c((nr_half %/% 2 + 1):nr_half, 1:(nr_half %/% 2)),
                     c((nc_half + 1):Nd, 1:nc_half)]
RDM <- Mod(sig_fft2)
RDM <- 10 * log10(RDM)

# Use the persp function (R's surf equivalent) to plot the output of 2DFFT
# and to show axis in both dimensions.
doppler_axis <- seq(-100, 100, length.out = Nd)
range_axis   <- seq(-200, 200, length.out = Nr / 2) * ((Nr / 2) / 400)

#open_device("02_rdm_surf")
#graphics::persp(doppler_axis, range_axis, t(RDM),
#                theta = 25, phi = 15, shade = 0.5, col = "lightblue",
#                border = NA, ticktype = "detailed",
#                xlab = "Velocity (m/s)", ylab = "Range (m)", zlab = "dB",
#                main = "Range Doppler Response")
#close_device()

## CFAR implementation

# Slide Window through the complete Range Doppler Map

# *TODO* :
# Select the number of Training Cells in both the dimensions.
# Tr <- ...
# Td <- ...

# *TODO* :
# Select the number of Guard Cells in both dimensions around the Cell under
# test (CUT) for accurate estimation.
# Gr <- ...
# Gd <- ...

# *TODO* :
# offset the threshold by SNR value in dB
# offset <- ...

# *TODO* :
# Create a vector to store noise_level for each iteration on training cells
noise_level <- 0

# *TODO* :
# design a loop such that it slides the CUT across range doppler map by
# giving margins at the edges for Training and Guard Cells.
# For every iteration sum the signal level within all the training
# cells. To sum convert the value from logarithmic to linear using
# 10^(x/10) (the db2pow equivalent). Average the summed values for all of
# the training cells used. After averaging convert it back to logarithm
# using 10*log10(...) (the pow2db equivalent). Further add the offset to
# it to determine the threshold. Next, compare the signal under CUT with
# this threshold. If the CUT level > threshold assign it a value of 1,
# else equate it to 0.

# Use RDM[i, j] as the matrix from the output of 2D FFT for implementing
# CFAR.

# *TODO* :
# The process above will generate a thresholded block, which is smaller
# than the Range Doppler Map as the CUT cannot be located at the edges of
# the matrix. Hence, few cells will not be thresholded. To keep the map
# size same set those values to 0.

# *TODO* :
# display the CFAR output using the persp function like we did for Range
# Doppler Response output.
#open_device("03_cfar_surf")
# graphics::persp(doppler_axis, range_axis, t(CFAR), ...)
#close_device()
