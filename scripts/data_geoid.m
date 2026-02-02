%% Process non-GHF datasets and save PNGs
close all; clear;

outdir = 'datasets_for_gmin';
addpath('data_ipics');
if ~exist(outdir,'dir'), mkdir(outdir); end

% --- PNG output folder ---
pngdir = fullfile(outdir, 'png');
if ~exist(pngdir, 'dir'), mkdir(pngdir); end

% %% === 0. Load Geoid ===
%load('geoid_N_egm2008.mat');  % lat_vec, lon_vec, N_grid

% Grid definition (in kilometers)
spacing_km = .1;
x_uniform = -200:spacing_km:800;   % 10001 points
y_uniform = -200:spacing_km:400;   % 6001 points
[Xgrid, Ygrid] = meshgrid(x_uniform * 1000, y_uniform * 1000);  % meters

% % Convert to lat/lon (for geoid sampling)
% projPS = projcrs(3031);
% [lat_target, lon_target] = projinv(projPS, Xgrid, Ygrid);
% lon_target(lon_target > 180) = lon_target(lon_target > 180) - 360;  % wrap
% 
% % Ensure geoid grid increases
% if any(diff(lon_vec) < 0)
%     [lon_vec, j] = sort(lon_vec); N_grid = N_grid(:, j);
% end
% if any(diff(lat_vec) < 0)
%     [lat_vec, i] = sort(lat_vec); N_grid = N_grid(i, :);
% end
% 
% % Interpolate geoid height (EGM2008)
% F_N = griddedInterpolant({lat_vec, lon_vec}, N_grid', 'linear', 'none');
% N = F_N(lat_target, lon_target);
% fprintf('Geoid height range: %.2f to %.2f m\n', min(N(:), [], 'omitnan'), max(N(:), [], 'omitnan'));
% 
% % Save geoid PNGs
% save_map_and_hist(N, x_uniform, y_uniform, ...
%     'Interpolated Geoid Height (EGM2008)', fullfile(pngdir,'geoid_N'), 'm');

% %% === 1. Ice Thickness ===
% H = interp_grd_to_target('data_ipics/icethk.xyz.grd','z','x','y',Xgrid,Ygrid,'m');
% assert_same_grid(H, Xgrid, Ygrid, 'H');
% 
% fprintf('Ice Thickness range: %.2f to %.2f m\n', min(H(:), [], 'omitnan'), max(H(:), [], 'omitnan'));
% save(fullfile(outdir,'coldex_icethk.mat'), 'H', 'Xgrid', 'Ygrid');
% 
% save_map_and_hist(H, x_uniform, y_uniform, ...
%     'Interpolated Ice Thickness', fullfile(pngdir,'icethk'), 'm');
% 
% %% === 2. Bed reflectivity ===
% BR = interp_grd_to_target('bedeco.xyz.grd','z','x','y',Xgrid,Ygrid,'km');
% assert_same_grid(BR, Xgrid, Ygrid, 'BR');
% 
% fprintf('Bed reflectivity stats (finite): min=%.3f, max=%.3f, NaNs=%d\n', ...
%     min(BR(:),[],'omitnan'), max(BR(:),[],'omitnan'), nnz(~isfinite(BR)));
% 
% save(fullfile(outdir,'coldex_bedeco.mat'), 'BR', 'Xgrid', 'Ygrid');
% 
% save_map_and_hist(BR, x_uniform, y_uniform, ...
%     'Interpolated Bed Reflectivity', fullfile(pngdir,'bedeco'), 'dB or a.u.');
% 
% % === 1c. Export combined x,y,icethk,br (EPSG:3031 meters) ===
% drop_nan_rows = true;
% 
% xv = Xgrid(:); yv = Ygrid(:);
% Hv = H(:);     BRv = BR(:);
% T = [xv, yv, Hv, BRv];
% 
% if drop_nan_rows
%     keep = isfinite(Hv) & isfinite(BRv);
%     T = T(keep, :);
% end
% 
% csv_out = fullfile(outdir,'xy_icethk_br.csv');
% fid = fopen(csv_out, 'w');
% fprintf(fid, "x_m,y_m,icethk_m,br\n");   % header
% fclose(fid);
% 
% % Correct argument order + append mode
% writematrix(single(T), csv_out, 'WriteMode','append');
% 
% mat_out = fullfile(outdir,'xy_icethk_br.mat');
% x_m = xv; y_m = yv; icethk = H; br = BR; 
% save(mat_out, 'x_m','y_m','icethk','br','-v7.3');
% 
% fprintf('Wrote combined table:\n  %s\nand MATLAB bundle:\n  %s\n', csv_out, mat_out);

% %% === 3. Surface Elevation (orthometric) ===
% S = interp_grd_to_target('data_ipics/srfelv.xyz.grd','z','x','y',Xgrid,Ygrid,'m');
% assert_same_grid(S, Xgrid, Ygrid, 'S (ellipsoidal)');
% S = S - N;   % ellipsoid → orthometric
% 
% 
% fprintf('Surface Elevation range: %.2f to %.2f m\n', min(S(:), [], 'omitnan'), max(S(:), [], 'omitnan'));
% save(fullfile(outdir,'coldex_srfelv.mat'), 'S', 'Xgrid', 'Ygrid');
% 
% % (Optional) set specific caxis before saving: e.g., caxis([prctile(S(:),1) prctile(S(:),99)])
% save_map_and_hist(S, x_uniform, y_uniform, ...
%     'Interpolated Surface Elevation (orthometric)', fullfile(pngdir,'srfelv'), 'm');

% === 3. Specularity Content == Prefer GeoTIFF; fallback to .xyz.grd
spec_tif = 'spec.xyz_val.tiff';      % <-- change filename if needed
spec_grd = 'spec.xyz_val.grd';

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

    xmin = min(xv); xmax = max(xv);
    ymin = min(yv); ymax = max(yv);
    inside = xs>=xmin & xs<=xmax & ys>=ymin & ys<=ymax;
    fprintf('[spec] target points inside source bbox: %d/%d (%.1f%%)\n', ...
        nnz(inside), numel(inside), 100*nnz(inside)/numel(inside));

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

figure; imagesc(x_uniform, y_uniform, Q);set(gca,'YDir','normal'); axis image; title('spec');

% Keep NaNs in MAT (do NOT use 9999 internally)
%Q(~isfinite(Q)) = NaN;

% Optional: create export-friendly 9999 version if you truly need it later
%Q_9999 = Q; 
%Q_9999(~isfinite(Q_9999)) = 9999;

save(fullfile(outdir,'specularity.mat'), 'Q', 'Xgrid', 'Ygrid','-v7.3');

save_map_and_hist(Q, x_uniform, y_uniform, ...
    'Specularity Content (regridded; NaNs→9999)', ...
    fullfile(pngdir,'specularity'), 'fraction', [0 0.5]);

% %% === 4. Bed Elevation (orthometric) ===
% B = interp_grd_to_target('data_ipics/bedelv.xyz.grd','z','x','y',Xgrid,Ygrid,'m');
% assert_same_grid(B, Xgrid, Ygrid, 'B (ellipsoidal)');
% B = B - N;
% 
% fprintf('Bed Elevation range: %.2f to %.2f m\n', min(B(:), [], 'omitnan'), max(B(:), [], 'omitnan'));
% save(fullfile(outdir,'coldex_bedelv.mat'), 'B', 'Xgrid', 'Ygrid');
% 
% save_map_and_hist(B, x_uniform, y_uniform, ...
%     'Interpolated Bed Elevation (orthometric)', fullfile(pngdir,'bedelv'), 'm');

% %% === 5. REMA Surface Elevation (orthometric) from 2 m GeoTIFF ===
% % Path to your clipped 2 m REMA tile (slightly larger than ROI)
% rema_tif = 'datasets_for_gmin/REMA_clipped.tif';   % <-- change to your actual filename
% 
% fprintf('\nReading 2 m REMA from %s\n', rema_tif);
% 
% % Read GeoTIFF
% try
%     [Z_rem_src, R_rem] = readgeoraster(rema_tif);   % Z_rem_src: MxN, R_rem: maprasterref / map reference
%     info_rem = georasterinfo(rema_tif);
% catch
%     % Older MATLAB fallback
%     [Z_rem_src, R_rem] = geotiffread(rema_tif);
%     info_rem = geotiffinfo(rema_tif);
% end
% 
% Z_rem_src = double(Z_rem_src);
% 
% % Handle NoData if available
% nodata_rem = [];
% try
%     if isfield(info_rem,'MissingDataIndicator')
%         nodata_rem = info_rem.MissingDataIndicator;
%     elseif isfield(info_rem,'NoData')
%         nodata_rem = info_rem.NoData;
%     end
% catch
% end
% if ~isempty(nodata_rem)
%     Z_rem_src(ismember(Z_rem_src, nodata_rem)) = NaN;
% end
% 
% % Build source x/y axes in source CRS
% [xv_src, yv_src] = worldGrid(R_rem);   % matrices
% x_src = xv_src(1, :);                  % 1 x Nx
% y_src = yv_src(:, 1);                  % Ny x 1
% 
% % Ensure ascending axes
% if numel(x_src) > 1 && x_src(2) < x_src(1)
%     x_src     = fliplr(x_src);
%     Z_rem_src = fliplr(Z_rem_src);
% end
% if numel(y_src) > 1 && y_src(2) < y_src(1)
%     y_src     = flipud(y_src);
%     Z_rem_src = flipud(Z_rem_src);
% end
% 
% % Set up CRS objects
% crs_tgt = projcrs(3031);   % your working EPSG:3031
% crs_src = [];
% try
%     if isfield(info_rem,'ProjectedCRS') && ~isempty(info_rem.ProjectedCRS)
%         crs_src = info_rem.ProjectedCRS;
%     elseif isfield(info_rem,'CoordinateReferenceSystem') && ~isempty(info_rem.CoordinateReferenceSystem)
%         crs_src = info_rem.CoordinateReferenceSystem;
%     end
% catch
% end
% if isempty(crs_src)
%     warning('REMA TIFF CRS not detected; assuming EPSG:3031.');
%     crs_src = projcrs(3031);
% end
% 
% % Interpolate REMA from 2 m grid onto your analysis grid
% F_rem = griddedInterpolant({y_src, x_src}, Z_rem_src, 'linear', 'none'); % note {row(y), col(x)}
% Z_rem_interp = F_rem(y_tgt_src, x_tgt_src);   % same shape as Xgrid,Ygrid
% 
% % Convert to orthometric by subtracting geoid (consistent with rest of script)
% REMA_srfelv = Z_rem_interp - N;
% 
% fprintf('REMA 2 m (regridded) elevation range: %.2f to %.2f m\n', ...
%     min(REMA_srfelv(:), [], 'omitnan'), max(REMA_srfelv(:), [], 'omitnan'));
% 
% save(fullfile(outdir,'REMA_srfelv.mat'), 'REMA_srfelv', 'Xgrid', 'Ygrid');
% 
% save_map_and_hist(REMA_srfelv, x_uniform, y_uniform, ...
%     'REMA Surface Elevation (2 m source → analysis grid)', ...
%     fullfile(pngdir,'REMA_srfelv'), 'm');

% %% === 6. Surface Temperature (.xyz in EPSG:3031 meters) ===
% Ts_raw = readmatrix('ts_yearly.csv.xyz', 'FileType','text');  % [X(m) Y(m) value]
% xts = Ts_raw(:,1); yts = Ts_raw(:,2); zts = Ts_raw(:,3) / 12; % monthly → annual mean
% [F_ts, kept_ts, frac_ts, nfin_ts] = make_interpolant(xts, yts, zts);
% Ts_interp = F_ts(Xgrid, Ygrid);
% 
% fprintf('Surface Temperature points: %d finite → %d unique (%.1f%% kept)\n', nfin_ts, kept_ts, 100*frac_ts);
% fprintf('Ts range: %.3f to %.3f (source units)\n', min(Ts_interp(:),[],'omitnan'), max(Ts_interp(:),[],'omitnan'));
% save(fullfile(outdir,'Ts_interp.mat'),'Ts_interp','Xgrid','Ygrid','-v7.3');
% 
% save_map_and_hist(Ts_interp, x_uniform, y_uniform, ...
%     'Interpolated Surface Temperature (annual mean)', fullfile(pngdir,'Ts_interp'), 'K (or °C)');
% 
% %% === 7. Mass Balance (.xyz in EPSG:3031 meters) — target: m ice eq yr^-1 ===
% % Input units: kg/m^2/y IE
% rho_ice = 917;     % kg m^-3  (bulk glacier ice)
% 
% Mb_raw = readmatrix('smbgl_yearly_ice.csv.xyz', 'FileType','text');  % [X(m) Y(m) value]
% xmb = Mb_raw(:,1); ymb = Mb_raw(:,2); zmb = Mb_raw(:,3) / rho_ice;
% 
% % build interpolant and grid on EPSG:3031 meters
% [F_mb, kept_mb, frac_mb, nfin_mb] = make_interpolant(xmb, ymb, zmb);
% Mb_interp = F_mb(Xgrid, Ygrid);  % now in m ice eq yr^-1
% 
% fprintf('[Mb] points: %d finite → %d unique (%.1f%% kept)\n', nfin_mb, kept_mb, 100*frac_mb);
% fprintf('[Mb] range (m ice eq yr^{-1}): %.4f .. %.4f\n', ...
%     min(Mb_interp(:),[],'omitnan'), max(Mb_interp(:),[],'omitnan'));
% 
% save(fullfile(outdir,'Mb_interp.mat'),'Mb_interp','Xgrid','Ygrid','-v7.3');
% 
% save_map_and_hist(Mb_interp, x_uniform, y_uniform, ...
%     'Interpolated Mass Balance (m ice eq yr^{-1})', fullfile(pngdir,'Mb_interp'), 'm ice eq yr^{-1}');
% 
% %% === Ice Velocity: antarctica_ice_velocity_450m_v2.nc ===
% vel_file = 'datasets_for_gmin/antarctica_ice_velocity_450m_v2.nc';
% 
% fprintf('\nReading ice velocity from %s\n', vel_file);
% 
% % read source grid (EPSG:3031 m)
% xv = ncread(vel_file, 'x');
% yv = ncread(vel_file, 'y');
% [Xv, Yv] = meshgrid(xv, yv);
% 
% VX0 = double(ncread(vel_file,'VX'));
% VY0 = double(ncread(vel_file,'VY'));
% 
% % Force [Ny x Nx]
% if isequal(size(VX0), [numel(xv) numel(yv)]), VX0 = VX0.'; end
% if isequal(size(VY0), [numel(xv) numel(yv)]), VY0 = VY0.'; end
% assert(isequal(size(VX0), [numel(yv) numel(xv)]), 'VX dims unexpected: %s', mat2str(size(VX0)));
% assert(isequal(size(VY0), [numel(yv) numel(xv)]), 'VY dims unexpected: %s', mat2str(size(VY0)));
% 
% % --- Ensure ascending axes ---
% if xv(2) < xv(1)
%     xv = flipud(xv);
%     VX0 = fliplr(VX0);
%     VY0 = fliplr(VY0);
% end
% if yv(2) < yv(1)
%     yv = flipud(yv);
%     VX0 = flipud(VX0);
%     VY0 = flipud(VY0);
% end
% 
% Fvx = griddedInterpolant({yv, xv}, VX0, 'linear', 'none');
% Fvy = griddedInterpolant({yv, xv}, VY0, 'linear', 'none');
% 
% vx = Fvx(Ygrid, Xgrid);
% vy = Fvy(Ygrid, Xgrid);
% speed = hypot(vx, vy);
% 
% assert_same_grid(speed, Xgrid, Ygrid, 'speed');
% assert_same_grid(vx, Xgrid, Ygrid, 'vx');
% assert_same_grid(vy, Xgrid, Ygrid, 'vy');
% 
% fprintf('Velocity speed range: %.2f to %.2f m/yr\n', ...
%     min(speed(:),[],'omitnan'), max(speed(:),[],'omitnan'));
% 
% save(fullfile(outdir,'mouginot_icevel.mat'), ...
%      'vx', 'vy', 'speed', 'Xgrid', 'Ygrid');
% 
% % maps and hist
% save_map_and_hist(speed, x_uniform, y_uniform, ...
%     'Ice velocity speed (m/yr)', fullfile(pngdir,'vel_speed'), 'm/yr');
% save_map_and_hist(vx, x_uniform, y_uniform, ...
%     'Ice velocity vx (m/yr)', fullfile(pngdir,'vel_vx'), 'm/yr');
% save_map_and_hist(vy, x_uniform, y_uniform, ...
%     'Ice velocity vy (m/yr)', fullfile(pngdir,'vel_vy'), 'm/yr');
% 
% assert_same_grid(N, Xgrid, Ygrid, 'Geoid N');
% assert_same_grid(H, Xgrid, Ygrid, 'Ice thickness H');
% assert_same_grid(BR, Xgrid, Ygrid, 'Bed reflectivity BR');
% assert_same_grid(S, Xgrid, Ygrid, 'Surface elevation S');
% assert_same_grid(Q, Xgrid, Ygrid, 'Specularity Q');
% assert_same_grid(B, Xgrid, Ygrid, 'Bed elevation B');
% assert_same_grid(Ts_interp, Xgrid, Ygrid, 'Surface temp Ts');
% assert_same_grid(Mb_interp, Xgrid, Ygrid, 'Mass balance Mb');
% assert_same_grid(speed, Xgrid, Ygrid, 'Velocity speed');
% 
% fprintf('All datasets saved in %s\nAll PNGs saved in %s\n', outdir, pngdir);

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

function Zt = interp_grd_to_target(ncfile, varZ, varX, varY, Xgrid, Ygrid, xyUnits)
% Interpolate a gridded NetCDF (x/y regular grid) onto target Xgrid/Ygrid (meters).
% xyUnits = 'km' if nc x,y are km; 'm' if meters.

    Z0 = double(ncread(ncfile, varZ));
    x  = double(ncread(ncfile, varX));
    y  = double(ncread(ncfile, varY));

    if strcmpi(xyUnits,'km')
        x = x(:)' * 1000;   % -> meters
        y = y(:)  * 1000;   % -> meters
    else
        x = x(:)';          % meters
        y = y(:);
    end

    % Force Z0 to [Ny x Nx] = [numel(y) x numel(x)]
    if isequal(size(Z0), [numel(x) numel(y)])
        Z0 = Z0.'; % swap to [y x]
    end
    if ~isequal(size(Z0), [numel(y) numel(x)])
        error('interp_grd_to_target: Z is %s; expected [%d %d] or [%d %d].', ...
              mat2str(size(Z0)), numel(y), numel(x), numel(x), numel(y));
    end

    % Ensure monotonic increasing axes for griddedInterpolant
    if numel(x)>1 && x(2)<x(1), x = fliplr(x); Z0 = fliplr(Z0); end
    if numel(y)>1 && y(2)<y(1), y = flipud(y); Z0 = flipud(Z0); end

    F = griddedInterpolant({y, x}, Z0, 'linear', 'none');  % {row(y), col(x)}
    Zt = F(Ygrid, Xgrid);
end

function assert_same_grid(Z, Xgrid, Ygrid, name)
    if ~isequal(size(Z), size(Xgrid)) || ~isequal(size(Z), size(Ygrid))
        error('%s size %s does not match target grid %s', name, mat2str(size(Z)), mat2str(size(Xgrid)));
    end
end

