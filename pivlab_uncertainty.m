% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% piv_uncertainty.m: Script to compute PIV uncertainty using PIVlab's built-in
% uncertainty functions: `uncertainty.disparity` and `uncertainty.error_window`.
% ------------------------------------------------------------------------------
% Experimental Fluid Mechanics Laboratory
% Department of Mechanical and Aerospace Engineering
% University of Central Florida, Orlando, FL, USA
% Author: Carlos Soto
% Edited: 2026-09-07
% ------------------------------------------------------------------------------
% MATLAB Requirements
% - Image Processing Toolbox.
% - PIVlab Add-On.
% - Parallel Processing Toolbox (optional but recommended).
% - Versions Tested: MATLAB 2025b with PIVlab 3.14.
% Required Files
% - Raw image file sequence (`.tif`).
% - Exported velocity field file sequence (`.txt`).
% - Exported PIVlab settings file (`.mat`).
%   - Alternatively, provide values manually in the script.
% Optional but Recommended Files
% - Exported PIVlab mask file (`.mat`). If masks are not available, the
%   following must be commented out from the script:
%   - In the "USER CONFIGURATION" section:
%     - The `file_msk` definition that specifies the mask file.
%   - In the "IMPORT SETTINGS AND MASKS" section:
%     - The `data_msk` definition that loads the masks.
%     - The `masks_all` definition that parses the masks.
%   - In the "UNCERTAINTY COMPUTATION" section:
%     - The masking code block in Step 4 of the `nr = 1:Nr` loop.
% - `ProgressBar.m`. This file defines the ProgressBar function used to track
%   the the timestep loop and should be placed in the same directory as this
%   script. It works for both single-threaded and multi-threaded execution.
%   If it is missing, the following must be commented out of the "UNCERTAINTY
%   COMPUTATION" section of the script:
%   - The `prog` definition immediately before the `nt = 1:Nt` loop.
%   - The `count(prog)` definition inside the `nt = 1:Nt` loop (final step
%     inside the loop).
% ------------------------------------------------------------------------------
% Note: File and directory paths should be configured to reflect the user's
% naming conventions and operating system.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure PIVlab's directory structure is added to the MATLAB path.
% To get the PIVlab toolbox installation directory run the following command:
% `disp(fullfile( ...
%   matlab.internal.addons.util.retrieveAddOnsInstallationFolder, ...
%   'Toolboxes', 'PIVlab'))`
addpath(genpath('/home/cesoto/MATLAB Add-Ons/Toolboxes/PIVlab'))

% Start a parallel pool for multi-threaded execution.
% If this is disabled, the timestep loop syntax in "UNCERTAINTY COMPUTATION"
% must be changed from `parfor = 1:Nt` to `for 1 = 1:Nt`.
current_pool = gcp('nocreate');
if isempty(current_pool)
  parpool('Processes')
else
  fprintf('\nParallel pool is ACTIVE with %d workers.\n', ...
    current_pool.NumWorkers);
end

% USER CONFIGURATION -----------------------------------------------------------

fprintf("\nDefining user settings...");

% Define PIV experiment parameters.
nc = 3;   % case to process
np = 1;   % plane to process
Nr = 5;   % number of runs
Nt = 180; % number of timesteps
fmt_t = "%0" + strlength(string(Nt)) + "d"; % format string for paths

% Define calibration error in pixels.
% This is based on human error in clicking the calibration grid points during
% processing. Typical values are between 0.5 px and 1 px.
err_cal = 0.5; % [px]

% Define t-penalty for 95% CI. For 5 runs (4 DOF), `tval_95=2.7764`.
tval_95 = 2.7764;

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

% Load PIVlab masks and settings.
data_msk = load(file_msk);
data_stg = load(file_stg);

% Parse masks.
masks_all = data_msk.masks_in_frame;

% Parse calibration settings.
% - If x_dir == 2, "x increases towards the left".
% - If y_dir == 2, "y increases towards the top".
dt        = str2double(data_stg.time_inp)/1000; % time step [s]
scale_xy  = data_stg.calxy; % calibration grid scaling factor [m/px]
scale_u   = data_stg.calu;  % u-velocity scaling factor [(m/s)/px]
scale_v   = data_stg.calv;  % v-velocity scaling factor [(m/s)/px]
L_grid    = str2double(data_stg.realdist)/1000; % Real calibration spacing [m]
N_px      = L_grid / scale_xy; % pixel calibration spacing [px]
x_dir     = data_stg.x_axis_direction; % x-direction
y_dir     = data_stg.y_axis_direction; % y-direction

% Parse cross-correlation settings.
win_size = round(str2double(data_stg.pass4val)); % final window size [px]
step_size = round(str2double(data_stg.step4)); % final step size [px]
img_thresh = str2double(data_stg.minintens); % image intensity threshold
img_size = data_stg.size_of_the_image; % image resolution [px]
img_y = img_size(1, 1); % image y-resolution [px]
img_x = img_size(1, 2); % image x-resolution [px]
Nx = floor(((img_x - win_size) / step_size)) + 1; % num x-points
Ny = floor(((img_y - win_size) / step_size)) + 1; % num y-points

% Parse image preprocessing settings.
clahe_enable  = data_stg.clahe_enable;
highp_enable  = data_stg.enable_highpass;
wiener_enable = data_stg.wienerwurst;
clahe_size    = str2double(data_stg.clahe_size);
highp_size    = str2double(data_stg.highp_size);
wiener_size   = str2double(data_stg.wienerwurstsize);

fprintf(' Done!\n');

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
  M_runs = zeros(Ny, Nx, Nr);
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

    % Get filepaths for image pair at this timestep/run.
    strA = sprintf(fmt_t, iA);
    strB = sprintf(fmt_t, iB);
    fold_img_run = compose(fold_img + "/img_c%d_p%d_r%d", nc, np, nr);
    file_imgA = compose(fold_img_run + "/img_c%d_p%d_r%d_%s.tif", ...
      nc, np, nr, strA);
    file_imgB = compose(fold_img_run + "/img_c%d_p%d_r%d_%s.tif", ...
      nc, np, nr, strB);

    % Get filepath for velocity field at this timestep/run.
    strT = sprintf(fmt_t, nt);
    fold_vel_run = compose(fold_vel + "/vel_c%d_p%d_r%d", nc, np, nr);
    file_vel = compose(fold_vel_run + "/vel_c%d_p%d_r%d_%s.txt", ...
      nc, np, nr, strT);

    % 2. Load raw image pair and preprocess ------------------------------------

    imgA_raw = double(imread(file_imgA));
    imgB_raw = double(imread(file_imgB));

    % HIGH-PASS FILTER (Subtract Gaussian blur to flatten background glow)
    if highp_enable == 1
      imgA_hp = imgA_raw - imfilter(imgA_raw, ...
        fspecial('gaussian', highp_size*4, highp_size), 'replicate');
      imgB_hp = imgB_raw - imfilter(imgB_raw, ...
        fspecial('gaussian', highp_size*4, highp_size), 'replicate');
    else
      imgA_hp = imgA_raw;
      imgB_hp = imgB_raw;
    end

    % WIENER FILTER (Denoise camera grain)
    if wiener_enable == 1
      imgA_w = wiener2(imgA_hp, [wiener_size, wiener_size]);
      imgB_w = wiener2(imgB_hp, [wiener_size, wiener_size]);
    else
      imgA_w = imgA_hp;
      imgB_w = imgB_hp;
    end

    % CLAHE FILTER (Normalize intensity and pop particle contrast)
    if clahe_enable == 1
      % Convert image temporarily back to 0-1 scale required by adapthisteq.
      imgA_norm = (imgA_w - min(imgA_w(:))) / (max(imgA_w(:)) - min(imgA_w(:)));
      imgB_norm = (imgB_w - min(imgB_w(:))) / (max(imgB_w(:)) - min(imgB_w(:)));
      % Apply CLAHE using your target tile size.
      % PIVlab sets NumTiles based on image dimension divided by clahe_size.
      tiles_y = round(img_y / clahe_size);
      tiles_x = round(img_x / clahe_size);
      imgA_prc = adapthisteq(imgA_norm, 'NumTiles', [tiles_y, tiles_x], ...
        'ClipLimit', 0.01) * 255;
      imgB_prc = adapthisteq(imgB_norm, 'NumTiles', [tiles_y, tiles_x], ...
        'ClipLimit', 0.01) * 255;
    else
      imgA_prc = imgA_w;
      imgB_prc = imgB_w;
    end

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
    M = hypot(U, V); % magnitude [m/s]

    % PIVlab exports columns sequentially. Reshape to matrix grid [Ny, Nx].
    X_grid = reshape(X, [Ny, Nx]);
    Y_grid = reshape(Y, [Ny, Nx]);
    U_grid = reshape(U, [Ny, Nx]);
    V_grid = reshape(V, [Ny, Nx]);
    M_grid = reshape(M, [Ny, Nx]);

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
    end
    if y_dir ~= 2
      X_grid = flipud(X_grid);
      Y_grid = flipud(Y_grid);
      U_grid = flipud(U_grid);
      V_grid = flipud(V_grid);
      M_grid = flipud(M_grid);
    end

    % 4. Compute displacement uncertainties ------------------------------------

    % Create window grid in pixels.
    first_center = win_size / 2; % first center is at half the window size
    grdx = ((0:(Nx-1)) * step_size) + first_center;
    grdy = (((0:(Ny-1)) * step_size) + first_center)';

    % Compute disparity components using PIVlab `uncertainty.disparity`
    % function.
    rsearch = 8;      % PIVlab default is 8
    weight = 'peaks'; % PIVlab default is 'peaks'
    [dispx, dispy, immult, peaks] = uncertainty.disparity( ...
      imgA_prc, imgB_prc, rsearch, weight);

    % Compute aggregate window uncertainty components using PIVlab
    % `uncertainty.error_window` function.
    ROI = [];             % leave empty to process entire frame
    nRmsLength = 1;       % PIVlab default is 1
    pos_weight = 'gauss'; % PIVlab default is 'gauss'
    warning('off', 'MATLAB:colon:nonIntegerIndex');
    [etotx, ebiasx, ermsx, ~, ~] = uncertainty.error_window( ...
      dispx, peaks, grdx, grdy, win_size, nRmsLength, pos_weight, ROI);
    [etoty, ebiasy, ermsy, ~, ~] = uncertainty.error_window( ...
      dispy, peaks, grdx, grdy, win_size, nRmsLength, pos_weight, ROI);
    warning('on', 'MATLAB:colon:nonIntegerIndex');

    % Get poly mask and convert to grid.
    [XG, YG] = meshgrid(grdx, grdy);
    mask_poly = masks_all{1, nt}{1, 2};
    mask_x = mask_poly(:, 1);
    mask_y = mask_poly(:, 2);
    mask_grid = inpolygon(XG, YG, mask_x, mask_y);

    % Apply mask to uncertainties.
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

    % Convert pixel displacement uncertainties to velocity.
    % Absolute value must be used because pixel displacements uncertainties
    % are unsigned scalar magnitudes.
    unc_disp_u_tr = etotx * abs(scale_u);
    unc_disp_v_tr = etoty * abs(scale_v);
    unc_disp_m_tr = hypot(unc_disp_u_tr, unc_disp_v_tr);

    % 5. Compute calibration and systematic uncertainties ----------------------

    % Compute calibration uncertainty.
    unc_cal_u_tr = U_grid * err_cal / N_px;
    unc_cal_v_tr = V_grid * err_cal / N_px;
    unc_cal_m_tr = M_grid * err_cal / N_px;

    % Compute systematic uncertainty.
    unc_sys_u_tr = hypot(unc_disp_u_tr, unc_cal_u_tr);
    unc_sys_v_tr = hypot(unc_disp_v_tr, unc_cal_v_tr);
    unc_sys_m_tr = hypot(unc_disp_m_tr, unc_cal_m_tr);

    % Store data at this timestep/run.
    X_runs(:, :, nr) = X_grid;
    Y_runs(:, :, nr) = Y_grid;
    U_runs(:, :, nr) = U_grid;
    V_runs(:, :, nr) = V_grid;
    M_runs(:, :, nr) = M_grid;
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

  % Compute max displacement uncertainty across runs.
  unc_disp_u_t = max(unc_disp_u_runs, [], 3, 'omitnan');
  unc_disp_v_t = max(unc_disp_v_runs, [], 3, 'omitnan');
  unc_disp_m_t = max(unc_disp_m_runs, [], 3, 'omitnan');

  % Compute max calibration uncertainty across runs.
  unc_cal_u_t = max(unc_cal_u_runs, [], 3, 'omitnan');
  unc_cal_v_t = max(unc_cal_v_runs, [], 3, 'omitnan');
  unc_cal_m_t = max(unc_cal_m_runs, [], 3, 'omitnan');

  % Compute max systematic uncertainty across runs.
  unc_sys_u_t = max(unc_sys_u_runs, [], 3, 'omitnan');
  unc_sys_v_t = max(unc_sys_v_runs, [], 3, 'omitnan');
  unc_sys_m_t = max(unc_sys_m_runs, [], 3, 'omitnan');

  % Compute mean grid and field data accross runs (for viz).
  X_avg_t = mean(X_runs, 3, 'omitnan');
  Y_avg_t = mean(Y_runs, 3, 'omitnan');
  U_avg_t = mean(U_runs, 3, 'omitnan');
  V_avg_t = mean(V_runs, 3, 'omitnan');
  M_avg_t = mean(M_runs, 3, 'omitnan');

  % Compute expanded random uncertainty across runs using 95% CI.
  U_std_t = std(U_runs, 0, 3, 'omitnan');
  V_std_t = std(V_runs, 0, 3, 'omitnan');
  M_std_t = std(M_runs, 0, 3, 'omitnan');
  unc_rand_u_t = tval_95 * U_std_t / sqrt(Nr);
  unc_rand_v_t = tval_95 * V_std_t / sqrt(Nr);
  unc_rand_m_t = tval_95 * M_std_t / sqrt(Nr);

  % Compute total combined systematic and expanded random uncertainty.
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

% Loop through timesteps and export uncertainties.
for nt = 1:Nt
  % Create filepath for current timestep.
  fxt = ".csv";
  strT = sprintf(fmt_t, nt);
  file_stat = compose(fold_stat + "/piv_stat_c%d_p%d_%s%s", nc, np, strT, fxt);

  % Reshape data at this time step to column vectors
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END: piv_uncertainty.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% CHECK: Plot uncertainty fields -----------------------------------------------

% Define uncertainty to plot
uncertainty_slice = unc_disp_m(:, :, 1);
% 1. Normalize your uncertainty grid between 0 and 1
min_val = min(uncertainty_slice(:));
max_val = max(uncertainty_slice(:));
normalized_grid = (uncertainty_slice - min_val) / (max_val - min_val);
% 2. Convert the normalized grid into an RGB image using a colormap matrix
cmap = hot(256);
rgb_image = ind2rgb(im2uint8(normalized_grid), cmap);
% 3. Use imwrite to dump the pixel data directly to a PNG
% This runs natively in a CLI with zero graphic requirements
imwrite(rgb_image, 'tmp/uncertainty_tmp.png');

% CHECK: All uncertainties should be positive ----------------------------------

% Define uncertainty to check
unc_check = unc_disp_v(:);
% Check that all values are positive
is_positive = (unc_check >= 0);
isnot_nan = ~isnan(uncertainty_slice);
all(is_positive(isnot_nan), 'all')

%}
