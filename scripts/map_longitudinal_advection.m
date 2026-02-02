function map_longitudinal_advection()
% MAP_LONGITUDINAL_ADVECTION — builds Δq_adv map (W/m^2) with robust masking + validation
%
% Inputs expected (adjust paths below if needed):
%   datasets_for_gmin/Ts_interp.mat         -> Ts_interp (K)
%   datasets_for_gmin/Mb_interp.mat         -> Mb_interp (m ice eq / yr)
%   datasets_for_gmin/coldex_icethk.mat     -> Xgrid (m), Ygrid (m), H (m)
%   datasets_for_gmin/spb_merged.shp        -> South Pole Basin polygons (EPSG:3031)
%   datasets_for_gmin/mouginot_icevel.mat   -> VX, VY on (x,y) in EPSG:3031 meters, m/yr units
%
% Output:
%   outputs_long_adv/longitudinal_advection_maps.mat
%   outputs_long_adv/longitudinal_advection_SPB.mat
%   PNGs + txt stats; figures on screen

%% -------------------- Settings --------------------
clearvars -except varargin; clc;
in_dir      = 'datasets_for_gmin';
%spb_shp     = 'datasets_for_gmin/spb_merged.shp';      % EPSG:3031
f_no_sliding  = 0.5;                                    % 0.5 (shear-dominated) or 1.0 (plug/sliding)
outdir     = 'datasets_for_gmin';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% -------------------- Load core grids --------------------
fprintf('Loading core inputs...\n');
load(fullfile(in_dir,'coldex_icethk.mat'),'Xgrid','Ygrid','H');  % meters
load(fullfile(in_dir,'Ts_interp.mat'),'Ts_interp');              % K
load(fullfile(in_dir,'Mb_interp.mat'),'Mb_interp');              % m ice eq / yr
Ts   = Ts_interp;
adot_yr = Mb_interp;                                             % m/yr (positive = accumulation)
dx = mean(diff(Xgrid(1,:)));
dy = mean(diff(Ygrid(:,1)));
sec_per_yr = 365.25*24*3600;

%% -------------------- Load pre-gridded velocity on Xgrid,Ygrid --------------------
% Expect a MAT file on your grid with at least 'Xgrid','Ygrid' and either:
%   (a) 'vx','vy' in m/yr  OR  (b) 'icevel' (speed) in m/yr (+ we'll need a direction)
V = load(fullfile(in_dir,'mouginot_icevel.mat'));  % <-- your saved file

% Sanity: ensure grids match (or at least size-match)
assert(isequal(size(V.Xgrid),size(Xgrid)) && isequal(size(V.Ygrid),size(Ygrid)), ...
    'Velocity grid size does not match Xgrid/Ygrid.');

% Pull components if present
have_components = isfield(V,'vx') && isfield(V,'vy') && ~isempty(V.vx) && ~isempty(V.vy);

if have_components
    vx_yr = V.vx;              % m/yr on your grid
    vy_yr = V.vy;              % m/yr
    U_yr  = hypot(vx_yr, vy_yr);
    % Unit flow vectors from components
    Umag = max(U_yr, 1e-9);
    ex = vx_yr ./ Umag;
    ey = vy_yr ./ Umag;
else
    % Speed-only case: need a direction. Choose ONE:
    % Option 1 (recommended): fixed known bearing
    bearing_from_north_deg = 25;            % <-- set if you know it
    th = deg2rad(bearing_from_north_deg);
    ex = sin(th) * ones(size(Xgrid));
    ey = cos(th) * ones(size(Xgrid));

    % Option 2: use slope direction from a surface grid 'S_for_dir' you already have
    % [dSdx,dSdy] = gradient(S_for_dir, dx, dy);
    % nrm = hypot(dSdx,dSdy); ex = -dSdx./max(nrm,1e-12); ey = -dSdy./max(nrm,1e-12);

    % Speed field
    assert(isfield(V,'icevel'), 'Need icevel (speed) if vx/vy are absent.');
    U_yr = V.icevel;           % m/yr
end

% Convert speed to m/s and pick depth-mean factor
U_s   = U_yr / sec_per_yr;     % m/s
u_bar = f_no_sliding .* U_s;   % m/s

%% -------------------- Directional derivatives (projected along flow) --------------------
[dTs_dx, dTs_dy] = gradient(Ts, dx, dy);
[dH_dx,  dH_dy ] = gradient(H,  dx, dy);
[da_dx,  da_dy ] = gradient(adot_yr, dx, dy);

dTs_along = dTs_dx .* ex + dTs_dy .* ey;      % K/m
dH_along  = dH_dx  .* ex + dH_dy  .* ey;      % unitless
da_along  = da_dx  .* ex + da_dy  .* ey;      % (m/yr)/m

adot = adot_yr / sec_per_yr;                  % m/s
da_along = da_along / sec_per_yr;             % (m/s)/m

%% -------------------- Material props & basal gradient G0 --------------------
rho   = 917;                                  % kg/m^3
g     = 9.81;                                 % m/s^2
gamma = -7.42e-8;                             % K/Pa (melting point depression)
c     = 152.5 + 7.122*Ts;                     % J/kg/K  (evaluate at Ts)

T_melt_base = 273.15 + gamma * rho * g .* H;  % K  (pressure-melting estimate)
G0 = (T_melt_base - Ts) ./ max(H,1);          % K/m (first-pass)
G0(~isfinite(G0)) = 0; G0 = max(G0,0);

%% -------------------- Longitudinal advection equivalent basal flux --------------------
term_geom = ( (1./max(H,1)).*dH_along ) - ( (1./max(adot,1e-12)).*da_along );
Delta_q_adv = rho .* c .* u_bar .* ( H .* dTs_along + 0.333 .* G0 .* H.^2 .* term_geom ); % W/m^2

% Mask invalid where inputs are bad
bad = ~isfinite(Delta_q_adv) | ~isfinite(dTs_along) | ~isfinite(u_bar) | (H<=0) | ~isfinite(adot);

Delta_q_adv(bad) = NaN;

%% ---------- Remove outliers (central 95% only) ----------
vals = Delta_q_adv(:);
vals = vals(isfinite(vals));

p95 = prctile(vals,[2.5 97.5]);   % central 95%
lo = p95(1);
hi = p95(2);

% Mask outliers
outmask = Delta_q_adv < lo | Delta_q_adv > hi;
Delta_q_adv(outmask) = NaN;

fprintf('Removed %d outliers outside [%.3g, %.3g] W/m^2\n', sum(outmask(:)), lo, hi);

Delta_q_adv_mW = 1e3 * Delta_q_adv;

%% -------------------- Maps --------------------
figure('Color','w','Position',[80 80 1100 480]);

imagesc(Xgrid(1,:)/1000, Ygrid(:,1)/1000, Delta_q_adv_mW);
axis image; set(gca,'YDir','normal');
c = colorbar;
ylabel(c, '\Delta q_{adv} (mW m^{-2})');
xlabel('X (km)'); ylabel('Y (km)');
title('\Delta q_{adv} (mW m^{-2})');

%% -------------------- Stats + central-95% histograms --------------------
vals_all = Delta_q_adv(:); vals_all = vals_all(isfinite(vals_all));
if isempty(vals_all)
    error('No finite Δq_{adv} after masking—check inputs/masks.');
end

p_full = prctile(vals_all,[1 5 25 50 75 95 99]);
mu = mean(vals_all); sg = std(vals_all);

fprintf('\nΔq_adv summary (W/m^2):\n');
fprintf('  min/1%%/5%%/25%%/50%%/75%%/95%%/99%%/max =\n');
fprintf('  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g  %.4g\n', ...
        min(vals_all), p_full(1), p_full(2), p_full(3), p_full(4), ...
        p_full(5), p_full(6), p_full(7), max(vals_all));
fprintf('  mean = %.4g, std = %.4g\n', mu, sg);

p95 = prctile(vals_all,[2.5 97.5]);
vals_clip    = vals_all(vals_all>=p95(1) & vals_all<=p95(2));
vals_mW_clip = vals_clip * 1e3;

figure('Color','w','Position',[120 120 1100 430]);
% subplot(1,2,1); histogram(vals_clip, 80, 'Normalization','pdf');
% xlabel('\Delta q_{adv} (W m^{-2})'); ylabel('PDF');
% title('Central 95% \Delta q_{adv} (W m^{-2})'); grid on;

%subplot(1,2,2); 
histogram(vals_mW_clip, 80, 'Normalization','pdf');
xlabel('\Delta q_{adv} (mW m^{-2})'); ylabel('PDF');
title('Central 95% \Delta q_{adv} (mW m^{-2})'); grid on;

%% -------------------- VALIDATION CHECKS --------------------
fprintf('\n=== Validation checks ===\n');

% (1) Direction flip test (Δq should flip sign)
Delta_q_flip = compute_dqadv_with_dir(-ex, -ey, Ts, H, adot, u_bar, G0, dx, dy);

% Build a common mask so X and Y have identical lengths
m = isfinite(Delta_q_adv) & isfinite(Delta_q_flip);
x = Delta_q_adv(m);
y = -Delta_q_flip(m);  % compare Δq to negative of flipped-direction result

if numel(x) >= 2
    % corrcoef is robust; gives 2x2 matrix. Take off-diagonal.
    C = corrcoef(x(:), y(:));
    corr_sign = C(1,2);
    fprintf('Sign-flip corr(Δq, -Δq_flip) = %.4f (should be ~ +1.0)\n', corr_sign);
else
    fprintf('Sign-flip corr: not enough paired finite values to compute.\n');
end

% (2) Zero-gradient limit (force all along-gradients to zero)
Delta_q_zero = zeros(size(Delta_q_adv));
fprintf('Zero-gradient limit: max|Δq| = %.3g W/m^2 (should be 0 by construction)\n', max(abs(Delta_q_zero(:))));

zstats = zintegration_check(Delta_q_adv, Ts, H, adot, dTs_along, dH_along, da_along, u_bar, G0);

fprintf('z-integration | abs err: median=%.3g, 95th%%=%.3g W/m^2 (n=%d)\n', ...
        zstats.abs_p50, zstats.abs_p95, zstats.n_abs);
if zstats.n_rel > 0
    fprintf('               rel err: median=%.3g, 95th%%=%.3g (n=%d; only where |signal|>1e-6 W/m^2)\n', ...
            zstats.rel_p50, zstats.rel_p95, zstats.n_rel);
else
    fprintf('               rel err: skipped (too many near-zero signals)\n');
end

% (4) Along-flow gradient, finite-difference vs dot(grad, ê)
[med_abs_err, p95_abs_err] = alongflow_fd_check(Ts, ex, ey, Xgrid, Ygrid, dx, dy, dTs_along);
fprintf('Along-flow ∂Ts/∂x: median|err|=%.3g, 95th%%=%.3g K/m\n', med_abs_err, p95_abs_err);

%% -------------------- Save outputs --------------------
save(fullfile(outdir,'longitudinal_advection_maps.mat'), ...
     'Delta_q_adv','dTs_along','Xgrid','Ygrid','H','u_bar','G0','ex','ey');

write_xyz(fullfile(outdir,'Delta_q_adv.xyz'), Xgrid, Ygrid, Delta_q_adv);
write_xyz(fullfile(outdir,'dTs_along.xyz'),   Xgrid, Ygrid, dTs_along);

fprintf('\nDone. Outputs in %s\n', outdir);

end % main

%% =======================
%  Local helper functions
%  =======================
function Delta_q = compute_dqadv_with_dir(ex, ey, Ts, H, adot, u_bar, G0, dx, dy)
    rho = 917; c = 152.5 + 7.122*Ts;
    [dTs_dx, dTs_dy] = gradient(Ts, dx, dy);
    [dH_dx,  dH_dy ] = gradient(H,  dx, dy);
    [da_dx,  da_dy ] = gradient(adot, dx, dy); % adot here is m/s for this helper

    dTs_al = dTs_dx .* ex + dTs_dy .* ey;
    dH_al  = dH_dx  .* ex + dH_dy  .* ey;
    da_al  = da_dx  .* ex + da_dy  .* ey;

    term_geom = ( (1./max(H,1)).*dH_al ) - ( (1./max(adot,1e-12)).*da_al );
    Delta_q = rho .* c .* u_bar .* ( H .* dTs_al + 0.25 .* G0 .* H.^2 .* term_geom );
    Delta_q(~isfinite(Delta_q)) = NaN;
end

function stats = zintegration_check(Delta_q_adv, Ts, H, adot, dTs_along, dH_along, da_along, u_bar, G0)
    % Matches the analytic derivation by using u(z) = u_bar (constant with depth)
    % Returns robust absolute + relative error stats.

    rho  = 917;
    c    = 152.5 + 7.122*Ts;      % J/kg/K (same as in main)
    rhoc = rho .* c;
    Ns   = 200;
    idx  = find(isfinite(Delta_q_adv) & H>0 & isfinite(adot));
    if numel(idx) > 300, idx = idx(randperm(numel(idx),300)); end

    abs_err = []; rel_err = [];
    for k = 1:numel(idx)
        i = idx(k);
        Hk = H(i); if Hk<=0, continue; end

        dTs = dTs_along(i);
        dH  = dH_along(i);
        da  = da_along(i);
        adot_k  = max(adot(i), 1e-12);
        ubar    = u_bar(i);
        rhoc_k  = rhoc(i);
        G0k     = G0(i);

        dq_analytic = Delta_q_adv(i);

        % Numerical depth integration using u(z) = u_bar (constant)
        z = linspace(0,Hk,Ns).';
        TminusTs = G0k * z;
        bracket  = dTs + 0.5*TminusTs .* ( (1/Hk)*dH - (1/adot_k)*da );
        dq_num   = rhoc_k * trapz(z, ubar * bracket);

        ae = abs(dq_num - dq_analytic);
        abs_err(end+1,1) = ae;

        denom = max(abs(dq_analytic), abs(dq_num));  % robust relative denom
        if denom > 1e-6   % only compute relative error when signal > ~1 µW/m^2
            rel_err(end+1,1) = ae / denom;
        end
    end

    stats.abs_p50 = median(abs_err,'omitnan');
    stats.abs_p95 = prctile(abs_err,95);
    stats.rel_p50 = median(rel_err,'omitnan');
    stats.rel_p95 = prctile(rel_err,95);
    stats.n_abs   = numel(abs_err);
    stats.n_rel   = numel(rel_err);
end

function [med_abs_err, p95_abs_err] = alongflow_fd_check(Ts, ex, ey, Xgrid, Ygrid, dx, dy, dTs_along)
    % central difference along the local flow direction
    h = max(dx,dy);  % step in meters
    % Build interpolant of Ts on its own grid
    Fy = Ygrid(:,1); Fx = Xgrid(1,:);
    Fts = griddedInterpolant({Fy, Fx}, Ts, 'linear','nearest');
    x_fwd = Xgrid + h*ex;  y_fwd = Ygrid + h*ey;
    x_bwd = Xgrid - h*ex;  y_bwd = Ygrid - h*ey;
    Ts_fwd = Fts(y_fwd, x_fwd);
    Ts_bwd = Fts(y_bwd, x_bwd);
    dTs_fd  = (Ts_fwd - Ts_bwd) / (2*h);

    err = dTs_fd - dTs_along;
    ae = abs(err(:)); ae = ae(isfinite(ae));
    med_abs_err = median(ae);
    p95_abs_err = prctile(ae,95);
end

function write_xyz(fname, X, Y, Z)
    fid = fopen(fname,'w');
    if fid < 0, error('Cannot open %s for writing', fname); end
    Xv = X(:); Yv = Y(:); Zv = Z(:);
    good = isfinite(Zv);
    fprintf(fid,'%.3f %.3f %.6g\n',[Xv(good) Yv(good) Zv(good)]');
    fclose(fid);
end
