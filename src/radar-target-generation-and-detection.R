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
target_range0 <- 80    # initial range (m)
target_vel    <- -70   # constant radial velocity (m/s), negative = approaching

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

## Signal generation and Moving Target simulation
# Running the radar scenario over the time.

for (i in seq_along(t)) {

  # *TODO* :
  # For each time stamp update the Range of the Target for constant velocity.
  # r_t[i] <- ...

  # *TODO* :
  # For each time sample we need update the transmitted and received signal.
  # Tx[i] <- ...
  # Rx[i] <- ...

  # *TODO* :
  # Now by mixing the Transmit and Receive generate the beat signal.
  # This is done by element wise matrix multiplication of Transmit and
  # Receiver Signal.
  # Mix[i] <- ...
}

## RANGE MEASUREMENT

# *TODO* :
# reshape the vector into Nr*Nd array. Nr and Nd here would also define the size
# of Range and Doppler FFT respectively.
# Mix2D <- matrix(Mix, nrow = Nr, ncol = Nd)

# *TODO* :
# run the FFT on the beat signal along the range bins dimension (Nr) and
# normalize.
# sig_fft <- fft(Mix2D[, 1]) / Nr

# *TODO* :
# Take the absolute value of FFT output
# sig_fft <- Mod(sig_fft)

# *TODO* :
# Output of FFT is double sided signal, but we are interested in only one side
# of the spectrum. Hence we throw out half of the samples.
# sig_fft <- sig_fft[1:(Nr / 2)]

# plotting the range
#open_device("01_range_fft")

# *TODO* :
# plot FFT output
# plot(..., type = "l", xlim = c(0, 200), ylim = c(0, 1),
#      xlab = "Range (m)", ylab = "Normalized amplitude")

#close_device()

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
