# Radar Target Generation and Detection — Project Report

Establish a complete radar target generation and detection pipeline: FMCW
waveform design, moving target simulation, range measurement (1st FFT),
Range-Doppler processing (2D FFT) and 2D CFAR-based target detection.

## Note on the implementation language (R instead of MATLAB)

This project is implemented in **R**, not MATLAB. Udacity's MATLAB license
for the Sensor Fusion course expired in **April 2026**, so the course
workspace no longer provides access to MATLAB. I tried to port the original MATLAB template (added here under `src/radar-target-generation-and-detection.m`) to R. The implementation uses base R only (`stats` + `graphics`) — no third-party packages are required.

## Project structure

| File | Purpose |
|---|---|
| `src/radar-target-generation-and-detection.m` | Original MATLAB template with open TODOs |
| `src/radar-target-generation-and-detection.R` | R port of the template.
Run the complete pipeline with:

```
Rscript src/radar-target-generation-and-detection.R
```

Plots are written as PNG files into `./plots`; in an interactive session
they are shown on screen.


## 1. FMCW Waveform Design

**Requirement:** Design the FMCW waveform from the given system requirements
and determine the bandwidth (B), the chirp time (Tchirp) and the chirp
slope. For the given requirements the calculated slope should be
approximately **2 × 10^13**.

System requirements:

- Carrier frequency: 77 GHz
- Maximum range: 200 m
- Range resolution: 1 m
- Maximum velocity: 100 m/s

Implementation: `src/radar-target-generation-and-detection.R`, section
*FMCW Waveform Generation* — the script computes and prints the parameters:

```
FMCW waveform design:
  Bandwidth  B      = 150.0 MHz
  Chirp time Tchirp = 7.333 us
  Slope      slope  = 2.0455e+13 Hz/s
```

Design:

- **Bandwidth** from the range resolution:
  `B = c / (2 * range_res) = 3e8 / 2 = 150 MHz`
- **Chirp time**, sized for the maximum range with a ~10% safety margin for
  the beat signal to be resolvable (factor 5.5 = 2 × 1.1 × 2.5, the common
  value used in the course):
  `Tchirp = 5.5 * (2 * max_range / c) = 5.5 * 1.33e-6 ≈ 7.33 µs`
- **Slope**:
  `slope = B / Tchirp = 150e6 / 7.33e-6 ≈ 2.045e13 Hz/s`

The slope of ≈ **2.045 × 10^13** meets the rubric requirement of approximately 2 × 10^13.

