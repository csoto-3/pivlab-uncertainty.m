---
documentclass: scrartcl
geometry: margin=1in
monofont: JuliaMono Nerd Font Mono
---

# PIVlab Uncertainty Computation Script

```text
Experimental Fluid Mechanics Laboratory
Department of Mechanical and Aerospace Engineering
University of Central Florida, Orlando, FL, USA
Author: Carlos Soto
Edited: 2026-09-10
```

This script (`pivlab_uncertainty.m`) was created to compute the PIV uncertainty
from an existing PIVlab dataset. In other words, if you already processed your
PIV data in PIVlab and need to compute the PIV uncertainty afterwards, this
script is for you.

DO NOT RELY ON THIS SCRIPT IF YOU HAVE NOT YET PROCESSED YOUR PIV DATA. JUST USE
THE UNCERTAINTY OPTIONS IN THE PIVLAB GUI DIRECTLY.

## Table of Contents

- [Requirements][hdr-req]
- [Uncertainty Calculation][hdr-unc]
- [Uncertainty Magnitude Computation][hdr-mag]

## Requirements

### MATLAB

- Image Processing Toolbox.
- PIVlab Add-On.
- Parallel Processing Toolbox (optional but recommended).
  - If parallel processing is disabled, all `parfor` loops must be changed to
    standard `for` loops.
- Versions Tested: MATLAB 2025b with PIVlab 3.14.

### Files

- Raw image file sequence (`.tif`).
- Exported velocity field file sequence (`.txt`).
- Exported PIVlab settings file (`.mat`).
  - Alternatively, provide values manually in the script.

> Note:
>
> Top-level directory paths should be configured to reflect the user's naming
> conventions and operating system.

#### Optional but Recommended Files

- Exported PIVlab mask file (`.mat`). If masks are not available, the
  following must be commented out from the script:
  - In the _USER CONFIGURATION_ section:
    - The `file_msk` definition that specifies the mask file.
  - In the _IMPORT SETTINGS AND MASKS_ section:
    - The `data_msk` definition that loads the masks.
    - The `masks_all` definition that parses the masks.
  - In the _UNCERTAINTY COMPUTATION_ section:
    - The masking code block in Step 4 of the `nr = 1:Nr` loop.
- `ProgressBar.m`. This file defines the `ProgressBar` function used to track
  the timestep loop and should be placed in the same directory as this script.
  It works for both single-threaded and multi-threaded execution.
  If it is missing, the following must be commented out of the script
  - All `prog` definitions placed before the start of `nt = 1:Nt` loops.
  - All `count(prog)` definitions inside the `nt = 1:Nt` loops (final step
    inside the loops).

## Uncertainty Calculation

$$
\vec{\sigma}_\mathrm{tot} =
\|\vec{\sigma}_\mathrm{rand} + \vec{\sigma}_\mathrm{sys}\|
$$

- $\vec{\sigma}_\mathrm{rand}=$ [random uncertainty][hdr-rand]
- $\vec{\sigma}_\mathrm{sys}=$ [systematic uncertainty][hdr-sys]
- Uncertainty magnitude is computed for by taking the
  [norm of the uncertainty subtype magnitudes][hdr-mag].

### Random Uncertainty

$$\vec{\sigma}_\mathrm{rand} = \vec{s}/\sqrt{N_r}$$

- $\vec{s}=$ standard deviation of velocity field across runs
- $N_r=$ number of runs
- Uncertainty magnitude is computed using the
  [error propagation method][hdr-mag].

### Systematic Uncertainty

$$
\vec{\sigma}_\mathrm{sys} =
\|\vec{\sigma}_\mathrm{disp} + \vec{\sigma}_\mathrm{cal}\|
$$

- $\vec{\sigma}_\mathrm{disp}=$ [displacement uncertainty][hdr-disp]
- $\vec{\sigma}_\mathrm{cal}=$ [calibration uncertainty][hdr-cal]
- Uncertainty components are computed for each run.
- Uncertainty magnitude is computed for each run by taking the
  [norm of the uncertainty subtype magnitudes][hdr-mag].
- The maximum of the 5 runs is then taken.

#### Displacement Uncertainty

$$\vec{\sigma}_\mathrm{disp}$$

- Uncertainty components are computed using PIVlab functions for each run.
- Uncertainty magnitude is computed for each run using the
  [error propagation method][hdr-mag].
- The maximum of the 5 runs is then taken.

Procedure

1. Image preprocessing (all optional)
   - Background Subtraction
     - Background image is precomputed using mean intensity of the image
       sequence. This is the default method in PIVlab.
   - CLAHE Filter
     - Filter settings taken from imported settings file or user-defined.
   - High-Pass Filter
     - Filter settings taken from imported settings file or user-defined.
   - Wiener Denoise Filter
     - Filter settings taken from imported settings file or user-defined.
   - Intensity Cap Filter -> NOT IMPLEMENTED
2. Disparity computation
   - Computed using the PIVlab `uncertainty.disparity` function.
   - Uses default PIVlab settings.
3. Pixel displacement uncertainty computation
   - Computed using the PIVlab `uncertainty.error_window` function.
   - Uses default PIVlab settings.
   - Takes disparity computation output as an input.
4. Application of pixel masks
   - Uncertainties in masked regions are set to `NaN`.
5. Uncertainty conversion from pixel displacement to velocity.
   - Scale factors taken from imported settings file or user-defined.
6. Filtering of interpolated velocity vectors.
   - Uncertainty of interpolated velocity vectors set to `NaN`.
   - Interpolated velocity vectors identified from `Vector Type [-]` column of
     imported velocity field data.

#### Calibration Uncertainty

$$
\vec{\sigma}_\mathrm{cal} =
\frac{\epsilon_\mathrm{cal}}{N_\mathrm{px}} \vec{u}
$$

- $\vec{u}=$ velocity vector
- $\epsilon_\mathrm{cal}=$ pixel calibration error
  - Result of human variance when clicking calibration grid points in the PIVlab
    GUI.
  - Typical values are between $0.5\,\mathrm{px}$ and $1.0\,\mathrm{px}$.
- $N_\mathrm{px}=$ number of pixels between calibration grid points
- Uncertainty components are computed for each run.
- Uncertainty magnitude is computed for each run using the
  [error propagation method][hdr-mag].
- The maximum of the 5 runs is then taken.

## Uncertainty Magnitude Computation

### Error Propagation Method

The error propagation method applies to directly-computed uncertainties:

- Displacement Uncertainty: $\sigma_{m,\mathrm{disp}}$
- Calibration Uncertainty: $\sigma_{m,\mathrm{cal}}$
- Random Uncertainty: $\sigma_{m,\mathrm{rand}}$

$$\sigma_m = \|\vec{u} + \vec{\sigma}\| / \|\vec{u}\|$$

- $\vec{u}=\langle u,v \rangle$ velocity vector
- $\vec{\sigma}=\langle \sigma_u,\sigma_v \rangle$ uncertainty vector

### Norm of Uncertainty Subtype Magnitudes

For uncertainties that are a combination of uncertainty subtypes, the error
propagation method does NOT apply. Instead, the norm of uncertainty subtype
magnitudes is computed:

Systematic Uncertainty Magnitude

$$
\sigma_{m,\mathrm{sys}} =
\|\sigma_{m,\mathrm{disp}} + \sigma_{m,\mathrm{cal}}\|
$$

Total Uncertainty Magnitude

$$
\sigma_{m,\mathrm{tot}} =
\|\sigma_{m,\mathrm{rand}} + \sigma_{m,\mathrm{sys}}\|
$$

<!-- prettier-ignore-start -->
[hdr-req]: #requirements
[hdr-unc]: #uncertainty-calculation
[hdr-rand]: #random-uncertainty
[hdr-sys]: #systematic-uncertainty
[hdr-disp]: #displacement-uncertainty
[hdr-cal]: #calibration-uncertainty
[hdr-mag]: #uncertainty-magnitude-computation
<!-- prettier-ignore-end -->
