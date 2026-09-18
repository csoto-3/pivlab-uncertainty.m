% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% piv_uncertainty.m: Script to compute PIV uncertainty using PIVlab's built-in
% uncertainty functions: `uncertainty.disparity` and `uncertainty.error_window`.
% ------------------------------------------------------------------------------
% Experimental Fluid Mechanics Laboratory
% Department of Mechanical and Aerospace Engineering
% University of Central Florida, Orlando, FL, USA
% Author: Carlos Soto
% Edited: 2026-09-18
% ------------------------------------------------------------------------------
% See README for documentation.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2026 Carlos Soto
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://gnu.org>.

% MATLAB CONFIGURATION ---------------------------------------------------------

% Make sure PIVlab's directory structure is added to the MATLAB path.
% To get the PIVlab toolbox installation directory run the following command:
% `disp(fullfile( ...
%   matlab.internal.addons.util.retrieveAddOnsInstallationFolder, ...
%   'Toolboxes', 'PIVlab'))`
addpath(genpath('/home/cesoto/MATLAB Add-Ons/Toolboxes/PIVlab'))

% Start a parallel pool for multi-threaded execution.
current_pool = gcp('nocreate');
if isempty(current_pool)
  parpool('Processes')
else
  fprintf('\nParallel pool is ACTIVE with %d workers.\n', ...
    current_pool.NumWorkers);
end

% USER CONFIGURATION -----------------------------------------------------------

% Define PIV experiment parameters.
% REQUIRED: (Nr, Nt). Used for computation loops.
% OPTIONAL: (nc, np, fmt_t). Used for path variables: (fold_, file_)
nc = 3;   % case to process
np = 8;   % plane to process
Nr = 5;   % number of runs
Nt = 180; % number of timesteps
fmt_t = "%0" + strlength(string(Nt)) + "d"; % format string for paths

fprintf("\nDefining user settings for case %d, plane %d...", nc, np);

% Define calibration error in pixels (ϵ_cal).
% This is based on human error in clicking the calibration grid points during
% processing. Typical values are between 0.5 px and 1 px.
err_cal = 0.5; % [px]

% Define optional features
subbkg_enable = true; % background subtraction
mask_enable = true;   % pixel masking (from PIVLAB .mat file)

% Define input paths.
script_dir = pwd();
fold_img = compose(script_dir + ...
  "/tmp/piv_prc/piv_prc_c%d_p%d/img_c%d_p%d", nc, np, nc, np);
fold_vel = compose(script_dir + ...
  "/tmp/piv_prc/piv_prc_c%d_p%d/vel_c%d_p%d", nc, np, nc, np);
file_stg = compose(script_dir + ...
  "/tmp/piv_prc/piv_prc_c%d_p%d/pivlab_c%d_p%d_stg.mat", nc, np, nc, np);
file_msk = compose(script_dir + ...
  "/tmp/piv_prc/piv_prc_c%d_p%d/pivlab_c%d_p%d_msk.mat", nc, np, nc, np);

% Define output directory and create it if it does not exist.
fold_stat = compose(script_dir + "/tmp/piv_stat/piv_stat_c%d_p%d", nc, np);
if ~isfolder(fold_stat)
  mkdir(fold_stat)
end

fprintf(' Done!\n');

% IMPORT SETTINGS AND MASKS ----------------------------------------------------

fprintf("\nImporting settings and masks...")

% Load and parse masks
% OPTIONAL: Comment out if mask_enable == false
data_msk = load(file_msk);
masks_all = data_msk.masks_in_frame;

% Load PIVlab settings.
% OPTIONAL: Comment out if settings will be defined manually.
data_stg = load(file_stg);

% Parse calibration settings (REQUIRED).
% - If x_dir == 2, "x increases towards the left".
% - If y_dir == 2, "y increases towards the top".
scale_xy  = data_stg.calxy; % calibration grid scaling factor [m/px]
scale_u   = data_stg.calu;  % u-velocity scaling factor [(m/s)/px]
scale_v   = data_stg.calv;  % v-velocity scaling factor [(m/s)/px]
L_grid    = str2double(data_stg.realdist)/1000; % Real calibration spacing [m]
N_px      = L_grid / scale_xy; % pixel calibration spacing [px] (Nₚₓ)
x_dir     = data_stg.x_axis_direction; % x-direction
y_dir     = data_stg.y_axis_direction; % y-direction

% Parse cross-correlation settings (REQUIRED)
win_size = round(str2double(data_stg.pass4val)); % final window size [px]
step_size = round(str2double(data_stg.step4)); % final step size [px]
img_thresh = str2double(data_stg.minintens); % image intensity threshold
img_size = data_stg.size_of_the_image; % image resolution [px]
img_y = img_size(1, 1); % image y-resolution [px]
img_x = img_size(1, 2); % image x-resolution [px]
Nx = floor(((img_x - win_size) / step_size)) + 1; % num x-points
Ny = floor(((img_y - win_size) / step_size)) + 1; % num y-points

% Parse image preprocessing settings.
% REQUIRED: `_enable` flags. Set to false to disable image preprocessing.
% OPTIONAL: `_size` flags. Will not be used if image processing is disabled.
clahe_enable  = data_stg.clahe_enable;
highp_enable  = data_stg.enable_highpass;
wiener_enable = data_stg.wienerwurst;
% intcap_enable = data_stg.enable_intenscap; % NOT IMPLEMENTED
clahe_size    = str2double(data_stg.clahe_size);
highp_size    = str2double(data_stg.highp_size);
wiener_size   = str2double(data_stg.wienerwurstsize);

fprintf(' Done!\n');

% GET FILE PATHS ---------------------------------------------------------------

% Get filepaths for images
items = dir(fold_img);
isDirFlags = [items.isdir] & ~ismember({items.name}, {'.', '..'});
subDirs = items(isDirFlags);
fullPaths = string(fullfile({subDirs.folder}, {subDirs.name}));
subdir_list = fullPaths(:);
file_img_all = strings(2 * Nt, Nr);
for nr = 1:Nr
  items = dir(subdir_list(nr));
  isFileFlags = ~ismember({items.name}, {'.', '..'});
  files = items(isFileFlags);
  fullPaths = string(fullfile({files.folder}, {files.name}));
  file_img_all(:, nr) = fullPaths;
end

% Get filepaths for velocity fields
items = dir(fold_vel);
isDirFlags = [items.isdir] & ~ismember({items.name}, {'.', '..'});
subDirs = items(isDirFlags);
fullPaths = string(fullfile({subDirs.folder}, {subDirs.name}));
subdir_list = fullPaths(:);
file_vel_all = strings(Nt, Nr);
for nr = 1:Nr
  items = dir(subdir_list(nr));
  isFileFlags = ~ismember({items.name}, {'.', '..'});
  files = items(isFileFlags);
  fullPaths = string(fullfile({files.folder}, {files.name}));
  file_vel_all(:, nr) = fullPaths;
end

% COMPUTE BACKGROUND IMAGES ----------------------------------------------------
% - This step should be performed if PIVlab background subtraction was used when
%   preprocessing the images.
% - Computes background image by computing mean intensity, which is the default
%   in PIVlab.
% - To check if they were computed correctly, write background images to file
%   with the following command, then open in system image viewer:
%   `imwrite(img_bkg(:, :, 1), 'tmp/uncertainty_tmp.png');`

img_bkg = zeros(img_y, img_x, Nr, 'uint8');
if subbkg_enable == 1
  fprintf('\nInitializing background image computation...\n');
  prog = ProgressBar(Nr, taskname='Computing background images:', ui='cli');
  for nr = 1:Nr
    [N_files, ~] = size(file_img_all);
    img_stack = zeros(img_y, img_x, N_files, 'uint8');
    parfor i = 1:N_files
      img_stack(:, :, i) = imread(file_img_all(i, nr));
    end
    img_bkg(:, :, nr) = uint8(mean(double(img_stack), 3));
    count(prog)
  end
end

% UNCERTAINTY COMPUTATION ------------------------------------------------------

fprintf('\nInitializing uncertainty computation...\n');

% Preallocate containers for timesteps.
X_avg = zeros(Ny, Nx, Nt);
Y_avg = zeros(Ny, Nx, Nt);
U_avg = zeros(Ny, Nx, Nt);
V_avg = zeros(Ny, Nx, Nt);
M_avg = zeros(Ny, Nx, Nt);
unc_cal_u   = zeros(Ny, Nx, Nt);
unc_cal_v   = zeros(Ny, Nx, Nt);
unc_cal_m   = zeros(Ny, Nx, Nt);
unc_disp_u  = zeros(Ny, Nx, Nt);
unc_disp_v  = zeros(Ny, Nx, Nt);
unc_disp_m  = zeros(Ny, Nx, Nt);
unc_sys_u   = zeros(Ny, Nx, Nt);
unc_sys_v   = zeros(Ny, Nx, Nt);
unc_sys_m   = zeros(Ny, Nx, Nt);
unc_rand_u  = zeros(Ny, Nx, Nt);
unc_rand_v  = zeros(Ny, Nx, Nt);
unc_rand_m  = zeros(Ny, Nx, Nt);
unc_tot_u   = zeros(Ny, Nx, Nt);
unc_tot_v   = zeros(Ny, Nx, Nt);
unc_tot_m   = zeros(Ny, Nx, Nt);

% Compute uncertainty across time steps and runs.
prog = ProgressBar(Nt, taskname='Computing uncertainties:', ui='cli');
parfor nt = 1:Nt

  % Compute A/B indices for image pairs.
  iA = 2 * nt - 1;
  iB = 2 * nt;

  % Preallocate containers for runs at this timestep.
  X_runs = zeros(Ny, Nx, Nr);
  Y_runs = zeros(Ny, Nx, Nr);
  U_runs = zeros(Ny, Nx, Nr);
  V_runs = zeros(Ny, Nx, Nr);
  unc_disp_u_runs = zeros(Ny, Nx, Nr);
  unc_disp_v_runs = zeros(Ny, Nx, Nr);
  unc_disp_m_runs = zeros(Ny, Nx, Nr);
  unc_cal_u_runs  = zeros(Ny, Nx, Nr);
  unc_cal_v_runs  = zeros(Ny, Nx, Nr);
  unc_cal_m_runs  = zeros(Ny, Nx, Nr);
  unc_sys_u_runs  = zeros(Ny, Nx, Nr);
  unc_sys_v_runs  = zeros(Ny, Nx, Nr);
  unc_sys_m_runs  = zeros(Ny, Nx, Nr);

  for nr = 1:Nr

    % 1. Get filepaths ---------------------------------------------------------

    file_imgA = file_img_all(iA, nr);
    file_imgB = file_img_all(iB, nr);
    file_vel = file_vel_all(nt, nr);

    % 2. Load raw image pair and apply preprocessing ---------------------------
    % Processing steps follow PIVlab defaults:
    %   A. Background subtraction
    %   B. CLAHE filter
    %   C. High-pass filter
    %   D. Wiener denoise filter
    %   E. Intensity capping

    % Load image pair and convert to `uint8` if necessary. (0-255)
    imgA_raw = imread(file_imgA);
    imgB_raw = imread(file_imgB);
    if ~isa(imgA_raw, 'uint8')
      imgA_raw = im2uint8(imgA_raw);
      imgB_raw = im2uint8(imgB_raw);
    end

    % Use raw image by default, these will be overwritten at each preprocessing
    % step that is enabled.
    imgA_prc = imgA_raw;
    imgB_prc = imgB_raw;

    % 2A. Background Subtraction
    if subbkg_enable == 1
      % Perform background subtraction.
      imgA_sub = imsubtract(imgA_prc, img_bkg(:, :, nr));
      imgB_sub = imsubtract(imgB_prc, img_bkg(:, :, nr));
      % Overwrite processed image pair.
      imgA_prc = imgA_sub;
      imgB_prc = imgB_sub;
    end

    % Convert to normalized double (0-1) for MATLAB preprocessing functions.
    imgA_prc = im2double(imgA_prc);
    imgB_prc = im2double(imgB_prc);

    % 2B. CLAHE Filter
    if clahe_enable == 1
      % Define CLAHE parameters.
      clip_limit = 0.01;        % PIVlab default is 0.01
      distribution = 'uniform'; % PIVlab default is 'uniform'
      % Compute dynamic grid tiles. PIVlab automatically calculates the
      % 'NumTiles' layout by dividing the total resolution by the target
      % window size.
      tiles_x = round(img_x / clahe_size);
      tiles_y = round(img_y / clahe_size);
      % Apply CLAHE parameters.
      imgA_clahe = adapthisteq(imgA_prc, 'NumTiles', [tiles_y, tiles_x], ...
        'ClipLimit', clip_limit, 'Distribution', distribution);
      imgB_clahe = adapthisteq(imgB_prc, 'NumTiles', [tiles_y, tiles_x], ...
        'ClipLimit', clip_limit, 'Distribution', distribution);
      % Overwrite processed image pair in regular 1-255 pixel space.
      imgA_prc = imgA_clahe;
      imgB_prc = imgB_clahe;
    end

    % 2C. High-Pass Filter
    if highp_enable == 1
      % Create Gaussian low-pass filter profile used by PIVlab.
      h_gaussian = fspecial('gaussian', highp_size, highp_size);
      % Generate blurred low-frequency baseline.
      imgA_blur = imfilter(imgA_prc, h_gaussian, 'replicate');
      imgB_blur = imfilter(imgB_prc, h_gaussian, 'replicate');
      % Subtract blur to execute high-pass filter.
      imgA_highp = imgA_prc - imgA_blur;
      imgB_highp = imgB_prc - imgB_blur;
      % Clamp output to zero to prevent negative intensity decimal drops.
      imgA_prc = max(imgA_highp, 0);
      imgB_prc = max(imgB_highp, 0);
    end

    % 2D. Wiener Denoise Filter
    if wiener_enable == 1
      imgA_wiener = wiener2(imgA_prc, [wiener_size, wiener_size]);
      imgB_wiener = wiener2(imgB_prc, [wiener_size, wiener_size]);
      imgA_prc = imgA_wiener;
      imgB_prc = imgB_wiener;
    end

    % 2E. Intensity Capping -> NOT IMPLEMENTED

    % Convert to 0-255 double for PIVlab uncertainty functions.
    imgA_prc = imgA_prc .* 255.0
    imgB_prc = imgB_prc .* 255.0

    % 3. Import PIVlab field data and reshape to grid --------------------------

    % Import field data, skipping the 3 default PIVlab header rows.
    opts = detectImportOptions(file_vel, 'NumHeaderLines', 3);
    opts.VariableNamingRule = 'preserve';
    txt_data = readmatrix(file_vel, opts);

    % Extract columns based on header layout.
    X = txt_data(:, 1); % x [m]
    Y = txt_data(:, 2); % y [m]
    U = txt_data(:, 3); % u [m/s]
    V = txt_data(:, 4); % v [m/s]
    M = hypot(U, V); % velocity magnitude [m/s]
    vtype = txt_data(:, 5); % vector type [-], to mask interpolated vectors

    % PIVlab exports columns sequentially. Reshape to matrix grid [Ny, Nx].
    X_grid = reshape(X, [Ny, Nx]);
    Y_grid = reshape(Y, [Ny, Nx]);
    U_grid = reshape(U, [Ny, Nx]);
    V_grid = reshape(V, [Ny, Nx]);
    M_grid = reshape(M, [Ny, Nx]);
    vtype_grid = reshape(vtype, [Ny, Nx]);

    % Flip the grid matrices based on PIVlab calibration.
    % - If x_dir == 2, "x increases towards the left" and the matrices must be
    %   flipped left-to-right.
    % - If y_dir ~= 2, "y increases towards the bottom" and the matrices must be
    %   flipped up-to-down.
    if x_dir == 2
      X_grid = fliplr(X_grid);
      Y_grid = fliplr(Y_grid);
      U_grid = fliplr(U_grid);
      V_grid = fliplr(V_grid);
      M_grid = fliplr(M_grid);
      vtype_grid = fliplr(vtype_grid);
    end
    if y_dir ~= 2
      X_grid = flipud(X_grid);
      Y_grid = flipud(Y_grid);
      U_grid = flipud(U_grid);
      V_grid = flipud(V_grid);
      M_grid = flipud(M_grid);
      vtype_grid = flipud(vtype_grid);
    end

    % 4. Compute displacement uncertainty (σ⃗_disp) -----------------------------

    % Create window grid in pixels.
    first_center = win_size / 2; % first center is at half the window size
    grdx = ((0:(Nx-1)) * step_size) + first_center;
    grdy = (((0:(Ny-1)) * step_size) + first_center)';

    % Compute disparity.
    rsearch = 8;      % PIVlab default is 8
    weight = 'peaks'; % PIVlab default is 'peaks'
    [dispx, dispy, immult, peaks] = uncertainty.disparity( ...
      imgA_prc, imgB_prc, rsearch, weight);

    % Compute pixel displacement uncertainty.
    ROI = [];             % leave empty to process entire frame
    nRmsLength = 1;       % PIVlab default is 1
    pos_weight = 'gauss'; % PIVlab default is 'gauss'
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    [etotx, ebiasx, ermsx, ~, ~] = uncertainty.error_window( ...
      dispx, peaks, grdx, grdy, win_size, nRmsLength, pos_weight, ROI);
    [etoty, ebiasy, ermsy, ~, ~] = uncertainty.error_window( ...
      dispy, peaks, grdx, grdy, win_size, nRmsLength, pos_weight, ROI);
    warning('on', 'MATLAB:colon:nonIntegerIndex');

    % Get pixel mask and convert to grid.
    [XG, YG] = meshgrid(grdx, grdy);
    mask_poly = masks_all{1, nt}{1, 2};
    mask_x = mask_poly(:, 1);
    mask_y = mask_poly(:, 2);
    mask_grid = inpolygon(XG, YG, mask_x, mask_y);

    % Apply pixel mask
    etotx(mask_grid)   = NaN;
    etoty(mask_grid)   = NaN;
    ebiasx(mask_grid)  = NaN;
    ebiasy(mask_grid)  = NaN;
    ermsx(mask_grid)   = NaN;
    ermsy(mask_grid)   = NaN;

    % Flip the pixel displacement matrices based on PIVlab calibration.
    % - If x_dir == 2, "x increases towards the left" and the matrices must be
    %   flipped left-to-right.
    % - If y_dir ~= 2, "y increases towards the bottom" and the matrices must be
    %   flipped up-to-down.
    if x_dir == 2
      etotx   = fliplr(etotx);
      etoty   = fliplr(etoty);
      ebiasx  = fliplr(ebiasx);
      ebiasy  = fliplr(ebiasy);
      ermsx   = fliplr(ermsx);
      ermsy   = fliplr(ermsy);
    end
    if y_dir ~= 2
      etotx   = flipud(etotx);
      etoty   = flipud(etoty);
      ebiasx  = flipud(ebiasx);
      ebiasy  = flipud(ebiasy);
      ermsx   = flipud(ermsx);
      ermsy   = flipud(ermsy);
    end

    % Apply vector filtering mask if enabled
    if mask_enable == 1
      mask_vtype = (vtype_grid == 2) | (vtype_grid == 0)
      etotx(mask_vtype)   = NaN;
      etoty(mask_vtype)   = NaN;
      ebiasx(mask_vtype)  = NaN;
      ebiasy(mask_vtype)  = NaN;
      ermsx(mask_vtype)   = NaN;
      ermsy(mask_vtype)   = NaN;
    end

    % Convert pixel displacement uncertainties to velocity (σ⃗_disp).
    unc_disp_u_tr = ebiasx * abs(scale_u);
    unc_disp_v_tr = ebiasy * abs(scale_v);
    unc_disp_m_tr = hypot(U_grid .* unc_disp_u_tr, V_grid .* unc_disp_v_tr) ...
      ./ M_grid;

    % 5. Compute calibration uncertainties (σ⃗_cal) -----------------------------
    unc_cal_u_tr = U_grid * err_cal / N_px;
    unc_cal_v_tr = V_grid * err_cal / N_px;
    unc_cal_m_tr = hypot(U_grid .* unc_cal_u_tr, V_grid .* unc_cal_v_tr) ...
      ./ M_grid;

    % 6. Compute systematic uncertainty (σ⃗_sys) --------------------------------
    unc_sys_u_tr = hypot(unc_disp_u_tr, unc_cal_u_tr);
    unc_sys_v_tr = hypot(unc_disp_v_tr, unc_cal_v_tr);
    unc_sys_m_tr = hypot(unc_disp_m_tr, unc_cal_m_tr);

    % 7. Store data at this timestep/run ---------------------------------------
    X_runs(:, :, nr) = X_grid;
    Y_runs(:, :, nr) = Y_grid;
    U_runs(:, :, nr) = U_grid;
    V_runs(:, :, nr) = V_grid;
    unc_disp_u_runs(:, :, nr) = unc_disp_u_tr;
    unc_disp_v_runs(:, :, nr) = unc_disp_v_tr;
    unc_disp_m_runs(:, :, nr) = unc_disp_m_tr;
    unc_cal_u_runs(:, :, nr)  = unc_cal_u_tr;
    unc_cal_v_runs(:, :, nr)  = unc_cal_v_tr;
    unc_cal_m_runs(:, :, nr)  = unc_cal_m_tr;
    unc_sys_u_runs(:, :, nr)  = unc_sys_u_tr;
    unc_sys_v_runs(:, :, nr)  = unc_sys_v_tr;
    unc_sys_m_runs(:, :, nr)  = unc_sys_m_tr;
  end

  % Compute mean displacement uncertainty across runs.
  unc_disp_u_t = mean(unc_disp_u_runs, 3, 'omitnan');
  unc_disp_v_t = mean(unc_disp_v_runs, 3, 'omitnan');
  unc_disp_m_t = mean(unc_disp_m_runs, 3, 'omitnan');

  % Compute mean calibration uncertainty across runs.
  unc_cal_u_t = mean(unc_cal_u_runs, 3, 'omitnan');
  unc_cal_v_t = mean(unc_cal_v_runs, 3, 'omitnan');
  unc_cal_m_t = mean(unc_cal_m_runs, 3, 'omitnan');

  % Compute mean systematic uncertainty across runs.
  unc_sys_u_t = mean(unc_sys_u_runs, 3, 'omitnan');
  unc_sys_v_t = mean(unc_sys_v_runs, 3, 'omitnan');
  unc_sys_m_t = mean(unc_sys_m_runs, 3, 'omitnan');

  % Compute mean field data across runs.
  X_avg_t = mean(X_runs, 3, 'omitnan');
  Y_avg_t = mean(Y_runs, 3, 'omitnan');
  U_avg_t = mean(U_runs, 3, 'omitnan');
  V_avg_t = mean(V_runs, 3, 'omitnan');
  M_avg_t = hypot(U_avg_t, V_avg_t);

  % Compute random uncertainty across runs (σ⃗_rand).
  U_std_t = std(U_runs, 0, 3, 'omitnan');
  V_std_t = std(V_runs, 0, 3, 'omitnan');
  unc_rand_u_t = U_std_t / sqrt(Nr);
  unc_rand_v_t = V_std_t / sqrt(Nr);
  unc_rand_m_t = hypot(U_avg_t .* unc_rand_u_t, V_avg_t .* unc_rand_v_t) ...
    ./ M_avg_t;

  % Compute total uncertainty (σ⃗_tot).
  unc_tot_u_t = hypot(unc_sys_u_t, unc_rand_u_t);
  unc_tot_v_t = hypot(unc_sys_v_t, unc_rand_v_t);
  unc_tot_m_t = hypot(unc_sys_m_t, unc_rand_m_t);

  % Store data at this timestep.
  X_avg(:, :, nt) = X_avg_t;
  Y_avg(:, :, nt) = Y_avg_t;
  U_avg(:, :, nt) = U_avg_t;
  V_avg(:, :, nt) = V_avg_t;
  M_avg(:, :, nt) = M_avg_t;
  unc_disp_u(:, :, nt)  = unc_disp_u_t;
  unc_disp_v(:, :, nt)  = unc_disp_v_t;
  unc_disp_m(:, :, nt)  = unc_disp_m_t;
  unc_cal_u(:, :, nt)   = unc_cal_u_t;
  unc_cal_v(:, :, nt)   = unc_cal_v_t;
  unc_cal_m(:, :, nt)   = unc_cal_m_t;
  unc_sys_u(:, :, nt)   = unc_sys_u_t;
  unc_sys_v(:, :, nt)   = unc_sys_v_t;
  unc_sys_m(:, :, nt)   = unc_sys_m_t;
  unc_rand_u(:, :, nt)  = unc_rand_u_t;
  unc_rand_v(:, :, nt)  = unc_rand_v_t;
  unc_rand_m(:, :, nt)  = unc_rand_m_t;
  unc_tot_u(:, :, nt)   = unc_tot_u_t;
  unc_tot_v(:, :, nt)   = unc_tot_v_t;
  unc_tot_m(:, :, nt)   = unc_tot_m_t;

  % Update progress bar.
  count(prog)
end

fprintf('\nDone!\n');

% EXPORT UNCERTAINTY FIELDS ---------------------------------------------------

fprintf("\nExporting uncertainty fields...");
for nt = 1:Nt

  % Define output file.
  fxt = ".csv";
  strT = sprintf(fmt_t, nt);
  file_stat = compose(fold_stat + "/piv_stat_c%d_p%d_%s%s", nc, np, strT, fxt);

  % Reshape data at this time step to column vectors.
  vec_x = reshape(X_avg(:, :, nt), [], 1);
  vec_y = reshape(Y_avg(:, :, nt), [], 1);
  vec_u = reshape(U_avg(:, :, nt), [], 1);
  vec_v = reshape(V_avg(:, :, nt), [], 1);
  vec_m = reshape(M_avg(:, :, nt), [], 1);
  vec_tot_u   = reshape(unc_tot_u(:, :, nt), [], 1);
  vec_tot_v   = reshape(unc_tot_v(:, :, nt), [], 1);
  vec_tot_m   = reshape(unc_tot_m(:, :, nt), [], 1);
  vec_sys_u   = reshape(unc_sys_u(:, :, nt), [], 1);
  vec_sys_v   = reshape(unc_sys_v(:, :, nt), [], 1);
  vec_sys_m   = reshape(unc_sys_m(:, :, nt), [], 1);
  vec_rand_u  = reshape(unc_rand_u(:, :, nt), [], 1);
  vec_rand_v  = reshape(unc_rand_v(:, :, nt), [], 1);
  vec_rand_m  = reshape(unc_rand_m(:, :, nt), [], 1);
  vec_cal_u   = reshape(unc_cal_u(:, :, nt), [], 1);
  vec_cal_v   = reshape(unc_cal_v(:, :, nt), [], 1);
  vec_cal_m   = reshape(unc_cal_m(:, :, nt), [], 1);
  vec_disp_u  = reshape(unc_disp_u(:, :, nt), [], 1);
  vec_disp_v  = reshape(unc_disp_v(:, :, nt), [], 1);
  vec_disp_m  = reshape(unc_disp_m(:, :, nt), [], 1);

  % Concatenate column vectors into matrices.
  mat_xy    = [vec_x, vec_y];
  mat_piv   = [vec_u, vec_v, vec_m];
  mat_tot   = [vec_tot_u, vec_tot_v, vec_tot_m];
  mat_sys   = [vec_sys_u, vec_sys_v, vec_sys_m];
  mat_rand  = [vec_rand_u, vec_rand_v, vec_rand_m];
  mat_cal   = [vec_cal_u, vec_cal_v, vec_cal_m];
  mat_disp  = [vec_disp_u, vec_disp_v, vec_disp_m];

  % Concatenate matrices into output table.
  mat_stat = [mat_xy, mat_piv, mat_tot, mat_sys, mat_rand, mat_cal, mat_disp];
  hdr_stat = {'X_avg', 'Y_avg', ...
    'U_avg', 'V_avg', 'M_avg', ...
    'U_unc_tot', 'V_unc_tot', 'M_unc_tot', ...
    'U_unc_sys', 'V_unc_sys', 'M_unc_sys', ...
    'U_unc_rand', 'V_unc_rand', 'M_unc_rand', ...
    'U_unc_cal', 'V_unc_cal', 'M_unc_cal', ...
    'U_unc_disp', 'V_unc_disp', 'M_unc_disp', ...
    };
  tbl_stat = array2table(mat_stat, 'VariableNames', hdr_stat);

  % Write table to file.
  writetable(tbl_stat, file_stat);
end

fprintf(' Done!\n');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: piv_uncertainty.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% CHECK: Plot uncertainty fields -----------------------------------------------

% Define uncertainty to plot
uncertainty_slice = unc_disp_m(:, :, 60);
% 1. Normalize your uncertainty grid between 0 and 1
min_val = min(uncertainty_slice(:));
max_val = max(uncertainty_slice(:));
normalized_grid = (uncertainty_slice - min_val) / (max_val - min_val);
% 2. Convert the normalized grid into an RGB image using a colormap matrix
cmap = hot(256);
rgb_image = ind2rgb(im2uint8(normalized_grid), cmap);
% 3. Use imwrite to dump the pixel data directly to a PNG
imwrite(rgb_image, 'tmp/uncertainty_tmp.png');

% CHECK: All uncertainties should be positive ----------------------------------

% Define uncertainty to check
uncertainty_slice = unc_disp_v(:);
% Check that all values are positive
is_positive = (uncertainty_slice >= 0);
isnot_nan = ~isnan(uncertainty_slice);
all(is_positive(isnot_nan), 'all')

%}
