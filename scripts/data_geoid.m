%% Process non-GHF datasets and save PNGs
close all; clear;

outdir = 'datasets_for_gmin';
addpath('data_ipics');
if ~exist(outdir,'dir'), mkdir(outdir); end

% --- PNG output folder ---
pngdir = fullfile(outdir, 'png');
if ~exist(pngdir, 'dir'), mkdir(pngdir); end

%% === 0. Load Geoid ===
load('geoid_N_egm2008.mat');  % lat_vec, lon_vec, N_grid

% Grid definition (in kilometers)
spacing_km = .1;
x_uniform = -200:spacing_km:1000;   % 1001 points
y_uniform = -200:spacing_km:400;   % 601 points
[Xgrid, Ygrid] = meshgrid(x_uniform * 1000, y_uniform * 1000);  % meters

% Convert to lat/lon (for geoid sampling)
projPS = projcrs(3031);
[lat_target, lon_target] = projinv(projPS, Xgrid, Ygrid);
lon_target(lon_target > 180) = lon_target(lon_target > 180) - 360;  % wrap

% Ensure geoid grid increases
if any(diff(lon_vec) < 0)
    [lon_vec, j] = sort(lon_vec); N_grid = N_grid(:, j);
end
if any(diff(lat_vec) < 0)
    [lat_vec, i] = sort(lat_vec); N_grid = N_grid(i, :);
end

% Interpolate geoid height (EGM2008)
F_N = griddedInterpolant({lat_vec, lon_vec}, N_grid', 'linear', 'none');
N = F_N(lat_target, lon_target);
fprintf('Geoid height range: %.2f to %.2f m\n', min(N(:), [], 'omitnan'), max(N(:), [], 'omitnan'));

% Save geoid PNGs
save_map_and_hist(N, x_uniform, y_uniform, ...
    'Interpolated Geoid Height (EGM2008)', fullfile(pngdir,'geoid_N'), 'm');

%% === 1. Ice Thickness ===
H_raw = ncread('data_ipics/icethk.xyz.grd', 'z')';
x_h   = ncread('data_ipics/icethk.xyz.grd', 'x')' / 1000;  % km
y_h   = ncread('data_ipics/icethk.xyz.grd', 'y')  / 1000;  % km
[Xh, Yh] = meshgrid(x_h, y_h);
Fh = scatteredInterpolant(Xh(:), Yh(:), H_raw(:), 'linear', 'none');
H  = Fh(Xgrid/1000, Ygrid/1000);  % Xgrid,Ygrid are meters → convert to km for Fh

fprintf('Ice Thickness range: %.2f to %.2f m\n', min(H(:), [], 'omitnan'), max(H(:), [], 'omitnan'));
save(fullfile(outdir,'coldex_icethk.mat'), 'H', 'Xgrid', 'Ygrid');

save_map_and_hist(H, x_uniform, y_uniform, ...
    'Interpolated Ice Thickness', fullfile(pngdir,'icethk'), 'm');

%% === 2. Bed reflectivity ===
BR_raw = ncread('bedeco.xyz.grd', 'z')';
x_br   = ncread('bedeco.xyz.grd', 'x')' / 1000;  % km
y_br   = ncread('bedeco.xyz.grd', 'y')  / 1000;  % km

% Ensure ascending axes just in case
if numel(x_br)>1 && x_br(2) < x_br(1)
    x_br  = fliplr(x_br); 
    BR_raw = fliplr(BR_raw);
end
if numel(y_br)>1 && y_br(2) < y_br(1)
    y_br  = fliplr(y_br); 
    BR_raw = flipud(BR_raw);
end

% Use griddedInterpolant (regular grid -> faster)
F_br = griddedInterpolant({y_br, x_br}, BR_raw, 'linear', 'none');
BR   = F_br(Ygrid/1000, Xgrid/1000);  % target grid in km for this interpolant

fprintf('Bed reflectivity stats (finite): min=%.3f, max=%.3f, NaNs=%d\n', ...
    min(BR(:),[],'omitnan'), max(BR(:),[],'omitnan'), nnz(~isfinite(BR)));

save(fullfile(outdir,'coldex_bedeco.mat'), 'BR', 'Xgrid', 'Ygrid');

save_map_and_hist(BR, x_uniform, y_uniform, ...
    'Interpolated Bed Reflectivity', fullfile(pngdir,'bedeco'), 'dB or a.u.');

% === 1c. Export combined x,y,icethk,br (EPSG:3031 meters) ===
drop_nan_rows = true;

xv = Xgrid(:); yv = Ygrid(:);
Hv = H(:);     BRv = BR(:);
T = [xv, yv, Hv, BRv];

if drop_nan_rows
    keep = isfinite(Hv) & isfinite(BRv);
    T = T(keep, :);
end

csv_out = fullfile(outdir,'xy_icethk_br.csv');
fid = fopen(csv_out, 'w');
fprintf(fid, "x_m,y_m,icethk_m,br\n");   % header
fclose(fid);

% Correct argument order + append mode
writematrix(single(T), csv_out, 'WriteMode','append');

mat_out = fullfile(outdir,'xy_icethk_br.mat');
x_m = xv; y_m = yv; icethk = H; br = BR; 
save(mat_out, 'x_m','y_m','icethk','br','-v7.3');

fprintf('Wrote combined table:\n  %s\nand MATLAB bundle:\n  %s\n', csv_out, mat_out);

%% === 3. Surface Elevation (orthometric) ===
S_raw = ncread('data_ipics/srfelv.xyz.grd', 'z')';
x_s   = ncread('data_ipics/srfelv.xyz.grd', 'x')' / 1000;
y_s   = ncread('data_ipics/srfelv.xyz.grd', 'y')  / 1000;
[Xs, Ys] = meshgrid(x_s, y_s);
Fs = scatteredInterpolant(Xs(:), Ys(:), S_raw(:), 'linear', 'none');
S  = Fs(Xgrid/1000, Ygrid/1000);
S  = S - N;   % ellipsoid → orthometric

fprintf('Surface Elevation range: %.2f to %.2f m\n', min(S(:), [], 'omitnan'), max(S(:), [], 'omitnan'));
save(fullfile(outdir,'coldex_srfelv.mat'), 'S', 'Xgrid', 'Ygrid');

% (Optional) set specific caxis before saving: e.g., caxis([prctile(S(:),1) prctile(S(:),99)])
save_map_and_hist(S, x_uniform, y_uniform, ...
    'Interpolated Surface Elevation (orthometric)', fullfile(pngdir,'srfelv'), 'm');

% === 3. Specularity Content == Prefer GeoTIFF; fallback to .xyz.grd
spec_tif = 'spec.xyz_val.tiff';      % <-- change filename if needed
spec_grd = 'data_ipics/spec.xyz.grd';

if exist(spec_tif,'file')
    % --- Read GeoTIFF ---
    try
        [Q_src, R_src] = readgeoraster(spec_tif);   % Q_src: MxN, R_src: maprasterref
        info = georasterinfo(spec_tif);
    catch
        % Older MATLAB: use geotiffread/info
        [Q_src, R_src] = geotiffread(spec_tif);
        info = geotiffinfo(spec_tif);
    end

    Q_src = double(Q_src);
    % Handle NoData if present
    nodata = [];
    try
        if isfield(info,'MissingDataIndicator')
            nodata = info.MissingDataIndicator;
        elseif isfield(info,'NoData')
            nodata = info.NoData;
        end
    catch, end
    if ~isempty(nodata)
        Q_src(ismember(Q_src, nodata)) = NaN;
    end

    % Source grid (native CRS)
    %[X_src, Y_src] = worldGrid(R_src);% not used directly; we build axes next
    [xv_src, yv_src] = worldGrid(R_src);   % matrices
    xv = xv_src(1, :);
    yv = yv_src(:, 1);

    % Ensure ascending axes
    if numel(xv)>1 && xv(2) < xv(1)
        xv = fliplr(xv);
        Q_src = fliplr(Q_src);
    end
    if numel(yv)>1 && yv(2) < yv(1)
        yv = flipud(yv);
        Q_src = flipud(Q_src);
    end

    % CRS objects
    crs_tgt = projcrs(3031);
    crs_src = [];
    try
        if isfield(info, 'ProjectedCRS') && ~isempty(info.ProjectedCRS)
            crs_src = info.ProjectedCRS;
        elseif isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
            crs_src = info.CoordinateReferenceSystem;
        end
    catch
        try
            if isfield(info,'PCS') && ~isempty(info.PCS)
                crs_src = projcrs(info.PCS);
            end
        catch
        end
    end
    if isempty(crs_src)
        warning('Specularity TIFF CRS not detected. Assuming EPSG:3031.');
        crs_src = projcrs(3031);
    end

    % Transform target grid to source CRS
    [lat_tgt, lon_tgt] = projinv(crs_tgt, Xgrid, Ygrid);
    [xs, ys] = projfwd(crs_src, lat_tgt, lon_tgt);

    % Interpolate (linear, no extrapolation)
    Fq = griddedInterpolant({yv, xv}, Q_src, 'linear', 'none'); % note {row(y), col(x)} order
    Q  = Fq(ys, xs);

else
    % --- Fallback: .xyz.grd path ---
    Q_raw = ncread(spec_grd, 'z')';
    x_q   = ncread(spec_grd, 'x')' / 1000;   % km
    y_q   = ncread(spec_grd, 'y')  / 1000;   % km
    [Xq, Yq] = meshgrid(x_q, y_q);
    Fq = scatteredInterpolant(Xq(:), Yq(:), double(Q_raw(:)), 'linear', 'none');
    Q  = Fq(Xgrid/1000, Ygrid/1000);
end

% Replace NaNs with 9999 (explicit convention)
Q(~isfinite(Q)) = 9999;

fprintf('Specularity range before fill: %.2f to %.2f\n', ...
    min(Q(Q~=9999), [], 'omitnan'), max(Q(Q~=9999), [], 'omitnan'));
fprintf('Specularity NaNs set to 9999: %d cells\n', nnz(Q==9999));

save(fullfile(outdir,'specularity.mat'), 'Q', 'Xgrid', 'Ygrid');

save_map_and_hist(Q, x_uniform, y_uniform, ...
    'Specularity Content (regridded; NaNs→9999)', ...
    fullfile(pngdir,'specularity'), 'fraction', [0 0.5]);

%% === 4. Bed Elevation (orthometric) ===
bedelv = ncread('data_ipics/bedelv.xyz.grd', 'z')';
x_b    = ncread('data_ipics/bedelv.xyz.grd', 'x')' / 1000;
y_b    = ncread('data_ipics/bedelv.xyz.grd', 'y')  / 1000;
[Xb, Yb] = meshgrid(x_b, y_b);
Fb = scatteredInterpolant(Xb(:), Yb(:), bedelv(:), 'linear', 'none');
B  = Fb(Xgrid/1000, Ygrid/1000);
B  = B - N;

fprintf('Bed Elevation range: %.2f to %.2f m\n', min(B(:), [], 'omitnan'), max(B(:), [], 'omitnan'));
save(fullfile(outdir,'coldex_bedelv.mat'), 'B', 'Xgrid', 'Ygrid');

save_map_and_hist(B, x_uniform, y_uniform, ...
    'Interpolated Bed Elevation (orthometric)', fullfile(pngdir,'bedelv'), 'm');

%% === 5. REMA Surface Elevation (orthometric) ===
load('data_ipics/REMA_200m_crop.mat');  % X_crop, Y_crop, Z_crop  (X/Y in km)
F_rem = scatteredInterpolant(X_crop(:) * 1000, Y_crop(:) * 1000, Z_crop(:), 'linear', 'none');
REMA_srfelv = F_rem(Xgrid, Ygrid) - N;

fprintf('REMA Elevation range: %.2f to %.2f m\n', min(REMA_srfelv(:), [], 'omitnan'), max(REMA_srfelv(:), [], 'omitnan'));
save(fullfile(outdir,'REMA_srfelv.mat'), 'REMA_srfelv', 'Xgrid', 'Ygrid');

save_map_and_hist(REMA_srfelv, x_uniform, y_uniform, ...
    'Interpolated REMA Surface Elevation (orthometric)', fullfile(pngdir,'REMA_srfelv'), 'm');

%% === 6. Surface Temperature (.xyz in EPSG:3031 meters) ===
Ts_raw = readmatrix('ts_yearly.csv.xyz', 'FileType','text');  % [X(m) Y(m) value]
xts = Ts_raw(:,1); yts = Ts_raw(:,2); zts = Ts_raw(:,3) / 12; % monthly → annual mean
[F_ts, kept_ts, frac_ts, nfin_ts] = make_interpolant(xts, yts, zts);
Ts_interp = F_ts(Xgrid, Ygrid);

fprintf('Surface Temperature points: %d finite → %d unique (%.1f%% kept)\n', nfin_ts, kept_ts, 100*frac_ts);
fprintf('Ts range: %.3f to %.3f (source units)\n', min(Ts_interp(:),[],'omitnan'), max(Ts_interp(:),[],'omitnan'));
save(fullfile(outdir,'Ts_interp.mat'),'Ts_interp');

save_map_and_hist(Ts_interp, x_uniform, y_uniform, ...
    'Interpolated Surface Temperature (annual mean)', fullfile(pngdir,'Ts_interp'), 'K (or °C)');

%% === 7. Mass Balance (.xyz in EPSG:3031 meters) — target: m ice eq yr^-1 ===
% Input units: kg/m^2/y IE
rho_ice = 917;     % kg m^-3  (bulk glacier ice)

Mb_raw = readmatrix('smbgl_yearly_ice.csv.xyz', 'FileType','text');  % [X(m) Y(m) value]
xmb = Mb_raw(:,1); ymb = Mb_raw(:,2); zmb = Mb_raw(:,3) / rho_ice;

% build interpolant and grid on EPSG:3031 meters
[F_mb, kept_mb, frac_mb, nfin_mb] = make_interpolant(xmb, ymb, zmb);
Mb_interp = F_mb(Xgrid, Ygrid);  % now in m ice eq yr^-1

fprintf('[Mb] points: %d finite → %d unique (%.1f%% kept)\n', nfin_mb, kept_mb, 100*frac_mb);
fprintf('[Mb] range (m ice eq yr^{-1}): %.4f .. %.4f\n', ...
    min(Mb_interp(:),[],'omitnan'), max(Mb_interp(:),[],'omitnan'));

save(fullfile(outdir,'Mb_interp.mat'),'Mb_interp');  % saved as m ice eq yr^-1

save_map_and_hist(Mb_interp, x_uniform, y_uniform, ...
    'Interpolated Mass Balance (m ice eq yr^{-1})', fullfile(pngdir,'Mb_interp'), 'm ice eq yr^{-1}');

%% === Ice Velocity: antarctica_ice_velocity_450m_v2.nc ===
vel_file = 'antarctica_ice_velocity_450m_v2.nc';

fprintf('\nReading ice velocity from %s\n', vel_file);

% read source grid (EPSG:3031 m)
xv = ncread(vel_file, 'x');
yv = ncread(vel_file, 'y');
[Xv, Yv] = meshgrid(xv, yv);

% Read velocities
VX = ncread(vel_file, 'VX')';
VY = ncread(vel_file, 'VY')';

% --- Ensure ascending axes ---
if xv(2) < xv(1)
    xv = flipud(xv);   % flip column vector
    VX = fliplr(VX);   % flip data horizontally
    VY = fliplr(VY);
end
if yv(2) < yv(1)
    yv = flipud(yv);
    VX = flipud(VX);
    VY = flipud(VY);
end

% speed magnitude
%icevel = hypot(VX, VY);

% interpolate to your target grid Xgrid,Ygrid (EPSG:3031 m)
Fvx = griddedInterpolant({yv, xv}, VX, 'linear', 'none');
Fvy = griddedInterpolant({yv, xv}, VY, 'linear', 'none');
vx = Fvx(Ygrid, Xgrid);
vy = Fvy(Ygrid, Xgrid);
speed = hypot(vx, vy);

fprintf('Velocity speed range: %.2f to %.2f m/yr\n', ...
    min(speed(:),[],'omitnan'), max(speed(:),[],'omitnan'));

save(fullfile(outdir,'mouginot_icevel.mat'), ...
     'vx', 'vy', 'speed', 'Xgrid', 'Ygrid');

% maps and hist
save_map_and_hist(speed, x_uniform, y_uniform, ...
    'Ice velocity speed (m/yr)', fullfile(pngdir,'vel_speed'), 'm/yr');
save_map_and_hist(vx, x_uniform, y_uniform, ...
    'Ice velocity vx (m/yr)', fullfile(pngdir,'vel_vx'), 'm/yr');
save_map_and_hist(vy, x_uniform, y_uniform, ...
    'Ice velocity vy (m/yr)', fullfile(pngdir,'vel_vy'), 'm/yr');


fprintf('All datasets saved in %s\nAll PNGs saved in %s\n', outdir, pngdir);

%% =======================
%  Local helper functions
%  (MATLAB allows local functions at end of script)
%  =======================

function [F, kept, frac, nfin] = make_interpolant(x, y, z)
    Mfin = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(Mfin); y = y(Mfin); z = z(Mfin);
    nfin = numel(x);
    if nfin==0
        warning('No finite points for interpolant.');
        F=[]; kept=0; frac=0; return;
    end
    [xyu, ia] = unique([x(:) y(:)], 'rows', 'stable');
    F = scatteredInterpolant(xyu(:,1), xyu(:,2), z(ia), 'linear', 'none');
    kept = numel(ia); frac = kept / nfin;
end

function save_map_and_hist(Z, xkm, ykm, titleStr, basePath, unitsStr, clims)
    % Map
    fig1 = figure('Color','w','Visible','off');
    imagesc(xkm, ykm, Z); axis xy image tight;
    if nargin >= 7 && ~isempty(clims), clim(clims); end
    cb = colorbar;
    if nargin>=6 && ~isempty(unitsStr), cb.Label.String = unitsStr; end
    xlabel('Easting (km)'); ylabel('Northing (km)');
    title(titleStr);
    savepng(fig1, [basePath '_map.png']);
    close(fig1);

    % Histogram (clip to clims if provided; this also drops 9999 in specularity)
    z = Z(isfinite(Z));
    if nargin >= 7 && ~isempty(clims)
        z = z(z >= clims(1) & z <= clims(2));
    end
    fig2 = figure('Color','w','Visible','off');
    histogram(z, 100);
    if nargin>=6 && ~isempty(unitsStr), xlabel(unitsStr); else, xlabel('value'); end
    ylabel('count');
    title(sprintf('%s — histogram (n=%d)', titleStr, numel(z)));
    savepng(fig2, [basePath '_hist.png']);
    close(fig2);
end


function savepng(figHandle, fname)
    try
        exportgraphics(figHandle, fname, 'Resolution', 220);
    catch
        set(figHandle, 'PaperPositionMode','auto');
        print(figHandle, fname, '-dpng', '-r220');
    end
end
