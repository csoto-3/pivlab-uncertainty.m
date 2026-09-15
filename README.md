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
Edited: 2026-09-11
```

This script (`pivlab_uncertainty.m`) was created to compute the PIV uncertainty
from an existing PIVlab dataset. In other words, if you already processed your
PIV data in PIVlab and need to compute the PIV uncertainty afterwards, this
script is for you.

DO NOT RELY ON THIS SCRIPT IF YOU HAVE NOT YET PROCESSED YOUR PIV DATA. JUST USE
THE UNCERTAINTY OPTIONS IN THE PIVLAB GUI DIRECTLY. See the _Sub-pixel
estimation, robustness, and uncertainty_ section in the [online PIVlab
manual][man-pivlab] for details.

## Table of Contents

- [Prerequisites][hdr-req]
- [Uncertainty Calculation][hdr-unc]
- [Error Propagation Method for Uncertainty Magnitude][hdr-erp]
- [RSS Method for Uncertainty Combinations][hdr-rss]

## Prerequisites

### MATLAB

- MATLAB (v2025b or newer).
- MATLAB Image Processing Toolbox.
- PIVlab Toolbox (v3.14 or newer) - Licensed under the MIT License. Available on
  the [MATLAB File Exchange][mat-pivlab] or [GitHub][gh-pivlab].
- MATLAB Parallel Processing Toolbox. _Optional but recommended._
  - If parallel processing is disabled, all `parfor` loops must be changed to
    standard `for` loops.
- `ProgressBar.m` - Licensed under the MIT License. Available on
  [GitHub][gh-prog].
  - This file defines the `ProgressBar` function used to track the time step
    loop and should be placed in the same directory as the
    `pivlab_uncertainty.m` script.
  - Works for both single-threaded and multi-threaded execution.

### Files

- Raw image file sequence (`.tif`).
- Exported velocity field file sequence (`.txt`).
- Exported PIVlab settings file (`.mat`).
  - Alternatively, provide values manually in the script.
- Exported PIVlab mask file (`.mat`). _Optional but recommended._
  - If masks are not available, the following must be commented out from the
    script before execution:
    - In the _USER CONFIGURATION_ section:
      - The `file_msk` definition that specifies the mask file.
    - In the _IMPORT SETTINGS AND MASKS_ section:
      - The `data_msk` definition that loads the masks.
      - The `masks_all` definition that parses the masks.
    - In the _UNCERTAINTY COMPUTATION_ section:
      - The masking code block in Step 4 of the `nr = 1:Nr` loop.
- _NOTE:_ Top-level directory paths should be configured in the
  `pivlab_uncertainty.m` script to reflect the user's naming conventions and
  operating system.

## Uncertainty Calculation

### Uncertainty Tree

- [Total Uncertainty][hdr-tot] $(\sigma_\mathrm{tot})$
  - [Random Uncertainty][hdr-rand] $(\sigma_\mathrm{rand})$
  - [Systematic Uncertainty][hdr-sys] $(\sigma_\mathrm{sys})$
    - [Calibration Uncertainty][hdr-cal] $(\sigma_\mathrm{cal})$
    - [Displacement Uncertainty][hdr-disp] $(\sigma_\mathrm{disp})$

### Total Uncertainty

$$
\vec{\sigma}_\mathrm{tot} =
\left\Vert\vec{\sigma}_\mathrm{rand} + \vec{\sigma}_\mathrm{sys}\right\Vert
$$

$$
\sigma_{M,\mathrm{tot}} =
\left\Vert\sigma_{M,\mathrm{rand}} + \sigma_{M,\mathrm{sys}}\right\Vert
$$

- $\vec{\sigma}_\mathrm{rand}$= [random uncertainty][hdr-rand]
- $\vec{\sigma}_\mathrm{sys}$= [systematic uncertainty][hdr-sys]
- Uncertainty magnitude is computed using the [RSS method][hdr-rss].

### Random Uncertainty

$$\vec{\sigma}_\mathrm{rand} = \frac{\vec{s}}{\sqrt{N_r}}$$

$$
\sigma_{M,\mathrm{rand}} =
\frac{\left\Vert\vec{u} \odot \vec{\sigma}_\mathrm{rand}\right\Vert}{\left\Vert\vec{u}\right\Vert}
$$

- $\vec{s}$= standard deviation of velocity field across runs
- $N_r$= number of runs
- Uncertainty magnitude is computed using the [error propagation
  method][hdr-erp].

### Systematic Uncertainty

$$
\vec{\sigma}_\mathrm{sys} =
\left\Vert\vec{\sigma}_\mathrm{cal} + \vec{\sigma}_\mathrm{disp}\right\Vert
$$

$$
\sigma_{M,\mathrm{sys}} =
\left\Vert\sigma_{M,\mathrm{cal}} + \sigma_{M,\mathrm{disp}}\right\Vert
$$

- $\vec{\sigma}_\mathrm{cal}$= [calibration uncertainty][hdr-cal]
- $\vec{\sigma}_\mathrm{disp}$= [displacement uncertainty][hdr-disp]
- Uncertainty components are computed for each run.
- Uncertainty magnitude is computed using the [RSS method][hdr-rss].
- The ensemble average across runs is then taken.

#### Calibration Uncertainty

$$
\vec{\sigma}_\mathrm{cal} =
\frac{\epsilon_\mathrm{cal}}{N_\mathrm{px}} \vec{u}
$$

$$
\sigma_{M,\mathrm{cal}} =
\frac{\left\Vert\vec{u} \odot \vec{\sigma}_\mathrm{cal}\right\Vert}{\left\Vert\vec{u}\right\Vert}
$$

- $\vec{u}$= velocity vector
- $\epsilon_\mathrm{cal}$= pixel calibration error
  - Result of human variance when selecting the calibration grid points PIVlab.
  - Typical values are between $0.5\,\mathrm{px}$ and $1.0\,\mathrm{px}$.
- $N_\mathrm{px}$= number of pixels between calibration grid points
- Uncertainty components are computed for each run.
- Uncertainty magnitude is computed for each run using the [error propagation
  method][hdr-erp].
- The ensemble average across runs is then taken.
- _NOTE_: This method assumes only two grid points were used for calibration. If
  multiple grid points were used, then the RMS residual error of the calibration
  equations must be used. The RMS method has not been implemented in the script.

#### Displacement Uncertainty

$$\vec{\sigma}_\mathrm{disp}$$

$$
\sigma_{M,\mathrm{disp}} =
\frac{\left\Vert\vec{u} \odot \vec{\sigma}_\mathrm{disp}\right\Vert}{\left\Vert\vec{u}\right\Vert}
$$

- Uncertainty components are computed using PIVlab functions for each run.
- Uncertainty magnitude is computed for each run using the [error propagation
  method][hdr-erp].
- The ensemble average across runs is then taken.

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
   - Intensity Cap Filter $\longrightarrow$ NOT IMPLEMENTED
2. Disparity computation
   - Computed using the PIVlab `uncertainty.disparity` function.
   - Uses default PIVlab settings.
3. Pixel displacement uncertainty computation
   - Computed using the PIVlab `uncertainty.error_window` function.
   - Uses default PIVlab settings.
   - Takes disparity computation output as an input.
   - The output `ebias` is used for all further uncertainty computations.
4. Application of pixel masks
   - Uncertainties in masked regions are set to `NaN`.
5. Uncertainty conversion from pixel displacement to velocity.
   - Scale factors taken from imported settings file or user-defined.
6. Filtering of interpolated velocity vectors.
   - Uncertainty of interpolated velocity vectors set to `NaN`.
   - Interpolated velocity vectors identified from `Vector Type [-]` column of
     imported velocity field data.

## Error Propagation Method for Uncertainty Magnitude

The error propagation method applies to uncertainties that are _directly
computed_. The equation shown below is a simplification of the standard error
propagation formula for a generic uncertainty magnitude.

$$
\sigma_M =
\frac{\sqrt{u^2\sigma_u^2 + v^2\sigma_v^2}}{\sqrt{u^2 + v^2}} =
\frac{\left\Vert\vec{u} \odot \vec{\sigma}\right\Vert}{\left\Vert\vec{u}\right\Vert}
$$

- $\vec{u}=\langle u,v \rangle$ velocity vector
- $\vec{\sigma}=\langle \sigma_u,\sigma_v \rangle$ uncertainty vector

### Derivation of the Simplified Error Propagation Formula

Let $\vec{u}=\langle u,v \rangle$ represent a single velocity vector and let
$\vec{\sigma}=\langle \sigma_u,\sigma_v \rangle$ represent a single
corresponding uncertainty vector. The velocity magnitude $M$ is then given by:

$$\left\Vert\vec{u}\right\Vert=\sqrt{u^2 + v^2}=M$$

The error propagation formula for the uncertainty magnitude $\sigma_M$ is then
given by:

$$
\sigma_M = \sqrt{
\left(\frac{\partial M}{\partial u}\sigma_u\right)^2 +
\left(\frac{\partial M}{\partial v}\sigma_v\right)^2}
$$

Recall: $M = \sqrt{u^2+v^2} = \left(u^2 + v^2\right)^{\frac{1}{2}}$

The partial derivatives $\partial M/\partial u$ and $\partial M/\partial v$ are
then computed as follows.

$$
\frac{\partial M}{\partial u} =
\frac{1}{2}\left(u^2 + v^2\right)^{-\frac{1}{2}} \cdot 2u =
\frac{u}{\sqrt{u^2 + v^2}} = \frac{u}{M}
$$

$$
\frac{\partial M}{\partial v} =
\frac{1}{2}\left(u^2 + v^2\right)^{-\frac{1}{2}} \cdot 2v =
\frac{v}{\sqrt{u^2 + v^2}} = \frac{v}{M}
$$

The partial derivatives are then substituted into the error propagation formula
as follows.

$$
\sigma_M = \sqrt{
\left(\frac{u}{M}\sigma_u\right)^2 + \left(\frac{v}{M}\sigma_v\right)^2 } =
\sqrt{\frac{1}{M^2}u^2\sigma_u^2 + \frac{1}{M^2}v^2\sigma_v^2} =
\frac{1}{M}\sqrt{u^2\sigma_u^2 + v^2\sigma_v^2}
$$

This can be further simplified into vector form using the _Hadamard Product_
operator $\odot$, resulting in the expression presented earlier in this section.

$$
\sigma_M = \frac{\sqrt{u^2\sigma_u^2 + v^2\sigma_v^2}}{\sqrt{u^2 + v^2}} =
\frac{\left\Vert\vec{u} \odot \vec{\sigma}\right\Vert}{\left\Vert\vec{u}\right\Vert}
$$

## Root Sum Square (RSS) Method for Uncertainty Combinations

For uncertainties that are a combination of precomputed uncertainties, the _RSS
method_ is used. The error has already been propagated in the precomputed
uncertainties. Reapplying error propagation would incorrectly inflate the
uncertainty. This method applies to combinations of precomputed vectors _and_
magnitudes.

$$
\sigma_\mathrm{combo} =
\sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_N^2} =
\left\Vert\sigma_1 + \sigma_2 + \cdots + \sigma_N\right\Vert
$$

## License

This project is licensed under the GNU General Public License v3.0 - see the
[LICENSE](LICENSE) file for details.

<!-- prettier-ignore-start -->
[mat-pivlab]: https://www.mathworks.com/matlabcentral/fileexchange/27659-pivlab-particle-image-velocimetry-piv-tool-with-gui
[gh-pivlab]: https://github.com/Shrediquette/PIVlab
[gh-prog]: https://github.com/elgar328/matlab-code-examples/tree/main/tools/ProgressBar
[man-pivlab]: https://www.pivlab.de/manual/pages/piv-settings.html#shared
[hdr-req]: #prerequisites
[hdr-unc]: #uncertainty-calculation
[hdr-tot]: #total-uncertainty
[hdr-rand]: #random-uncertainty
[hdr-sys]: #systematic-uncertainty
[hdr-disp]: #displacement-uncertainty
[hdr-cal]: #calibration-uncertainty
[hdr-erp]: #error-propagation-method-for-uncertainty-magnitude
[hdr-rss]: #root-sum-square-rss-method-for-uncertainty-combinations
<!-- prettier-ignore-end -->
