%%% Process GHF Models (auto-CRS, unit-aware, orientation-safe) %%%

clearvars; close all;

%% ------------------ Toggles ------------------
doModels = struct( ...
    'Stal',     true, ...
    'Hazzard',  true, ...
    'An',       true, ...
    'FoxMaule', true, ...
    'Martos',   true, ...
    'Shen',     true, ...
    'Losing',   true, ...
    'save',     true ...
);

%% ------------------ Target grid (EPSG:3031) ------------------
tgrid = load('datasets_for_gmin/coldex_icethk.mat','Xgrid','Ygrid');
Xgrid = tgrid.Xgrid;
Ygrid = tgrid.Ygrid;     % meters, EPSG:3031
proj  = projcrs(3031);

qa_outdir = fullfile('datasets_for_gmin','qa_maps');
if ~exist(qa_outdir,'dir'), mkdir(qa_outdir); end

%% ------------------ Helpers ------------------
function Zs = ensure_grid_shape(Z, xv, yv)
    % Try to coerce Z to [numel(yv) x numel(xv)], otherwise return Z unchanged.
    ny = numel(yv); nx = numel(xv);
    [zy, zx] = size(Z);

    if zy == ny && zx == nx
        Zs = Z; return;
    end
    if zy == nx && zx == ny
        Zs = Z.'; return;  % pure transpose solves it
    end

    % If neither matches exactly, clip to the common extents safely.
    % This avoids index-out-of-bounds when some files have padding rows/cols.
    cy = min(zy, ny); cx = min(zx, nx);
    warning('ensure_grid_shape: size mismatch (Z=%dx%d vs axes=%dx%d). Clipping to %dx%d.', ...
        zy, zx, ny, nx, cy, cx);
    Zs = Z(1:cy, 1:cx);
end

function z = to_mW_if_W(filename, varname, z)
    u = '';
    try u = ncreadatt(filename, varname, 'units'); catch, end
    s = lower(string(u));
    if contains(s,'w') && ~contains(s,'mw')
        z = z * 1000;  % W/m^2 -> mW/m^2
        fprintf('[units] %s:%s converted W→mW\n', filename, varname);
    else
        fprintf('[units] %s:%s assumed mW\n', filename, varname);
    end
end

function A = clean_grid(A, file, varname)
    fv = []; mv = [];
    try fv = ncreadatt(file, varname, '_FillValue');   catch, end
    try mv = ncreadatt(file, varname, 'missing_value'); catch, end
    bad = [];
    if isnumeric(fv), bad = [bad; fv(:)]; end
    if isnumeric(mv), bad = [bad; mv(:)]; end
    bad = bad(isfinite(bad));
    M = ~isfinite(A);
    if ~isempty(bad), M = M | ismember(A, bad); end
    A(M) = NaN;
end

% ---------- AUTO CRS / UNITS DETECTION ----------
function [mode, ux, uy] = detect_axis_mode(file, xname, yname, xv, yv)
% mode ∈ {'degrees','meters','kilometers','unknown'}
    %mode = 'unknown'; 
    ux=''; uy='';
    try ux = lower(string(ncreadatt(file, xname, 'units'))); catch, end
    try uy = lower(string(ncreadatt(file, yname, 'units'))); catch, end
    % CF-Convention hints
    sx=''; sy='';
    try sx = lower(string(ncreadatt(file, xname, 'standard_name'))); catch, end
    try sy = lower(string(ncreadatt(file, yname, 'standard_name'))); catch, end
    xn = lower(string(xname)); yn = lower(string(yname));

    % 1) Attribute / name based
    if contains(ux,'degree') || contains(uy,'degree') || ...
       contains(sx,'longitude') || contains(sy,'latitude') || ...
       ismember(xn,["lon","longitude"]) || ismember(yn,["lat","latitude"])
        mode = 'degrees'; return;
    end
    if contains(ux,'kilometer') || contains(ux,'km') || contains(uy,'kilometer') || contains(uy,'km')
        mode = 'kilometers'; return;
    end
    if contains(ux,'metre') || contains(ux,'meter') || contains(uy,'metre') || contains(uy,'meter')
        mode = 'meters'; return;
    end

    % 2) Numeric heuristics (robust)
    ax = max(abs(xv)); ay = max(abs(yv));
    if ax <= 370 && ay <= 95
        mode = 'degrees'; return;
    end
    if ax > 5e3 || ay > 5e3
        mode = 'meters'; return;
    end

    % 3) Spacing heuristic
    dx = median(abs(diff(xv))); dy = median(abs(diff(yv)));
    if dx < 5 && dy < 5 && ax <= 720 && ay <= 180
        mode = 'degrees';
    else
        mode = 'meters';
    end
end

function [xv_m, yv_m, decided] = normalize_axes_units(file, xname, yname, xv, yv)
    [mode, ux, uy] = detect_axis_mode(file, xname, yname, xv, yv);
    decided = mode;
    switch mode
        case 'degrees'
            xv_m = xv; yv_m = yv;  % degrees for now (project later)
        case 'kilometers'
            xv_m = xv * 1000; yv_m = yv * 1000;
            fprintf('[axes] %s:(%s,%s) km→m\n', file, xname, yname);
        otherwise % meters or unknown→assume meters
            xv_m = xv; yv_m = yv;
    end
    % Log decision (only once)
    fprintf('[axes] %s:(%s,%s) detected=%s (ux="%s", uy="%s")\n', ...
        file, xname, yname, string(decided), string(ux), string(uy));
end

function [Xps, Yps, decided] = to3031_auto(file, xname, yname, xv, yv, proj)
    [xv_n, yv_n, decided] = normalize_axes_units(file, xname, yname, xv, yv);
    if strcmpi(decided,'degrees')
        [Xps, Yps] = to3031(xv_n, yv_n, true, proj);
    else
        [Xps, Yps] = to3031(xv_n, yv_n, false, proj);
    end
end
% ------------------------------------------------

function [xps,yps] = to3031(xv, yv, assume_deg, proj)
    if assume_deg
        [lon, lat] = meshgrid(xv, yv);
        [xps, yps] = projfwd(proj, lat, lon);
    else
        [xps, yps] = meshgrid(xv, yv);
    end
end

function [F, kept, frac, nfin] = make_interpolant(x, y, z)
    Mfin = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(Mfin); y = y(Mfin); z = z(Mfin);
    nfin = numel(x);
    if nfin==0
        warning('No finite points for interpolant. Returning empty.');
        F = []; kept = 0; frac = 0; return;
    end
    [xyu, ia] = unique([x(:) y(:)], 'rows', 'stable');
    F = scatteredInterpolant(xyu(:,1), xyu(:,2), z(ia), 'linear', 'none');
    kept = numel(ia); frac = kept / nfin;
end

function qa_map(Z, Xgrid, Ygrid, titleStr, out_png)
    figure('Color','w');
    imagesc([min(Xgrid(1,:)) max(Xgrid(1,:))], [min(Ygrid(:,1)) max(Ygrid(:,1))], Z);
    set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); colorbar;
    title(sprintf('%s\n[min=%.3g max=%.3g finite=%d/%d]', titleStr, ...
          min(Z(:),[],'omitnan'), max(Z(:),[],'omitnan'), nnz(isfinite(Z)), numel(Z)));
    xlabel('X (m, EPSG:3031)'); ylabel('Y (m, EPSG:3031)');
    try exportgraphics(gcf, out_png, 'Resolution', 180); catch, end
    close all;
end

function qa_hist(Z, titleStr, out_png)
    z = Z(isfinite(Z));
    figure('Color','w'); histogram(z, 100);
    title(sprintf('%s — histogram (n=%d)', titleStr, numel(z)));
    xlabel('mW m^{-2}'); ylabel('count');
    try exportgraphics(gcf, out_png, 'Resolution', 180); catch, end
    close all;
end

function coverage_check(name, xps, yps, Xgrid, Ygrid)
    in = xps(:)>=min(Xgrid(1,:)) & xps(:)<=max(Xgrid(1,:)) & ...
         yps(:)>=min(Ygrid(:,1)) & yps(:)<=max(Ygrid(:,1));
    fprintf('[cover] %s: points inside target bbox = %d/%d (%.1f%%)\n', ...
        name, nnz(in), numel(in), 100*nnz(in)/numel(in));
end

% Source-map helpers (same as previous answer; omitted here for brevity)
function qa_map_rect_source(xv, yv, Z, assume_deg, proj, titleStr, out_png)
    % Coerce size
    Zc = ensure_grid_shape(Z, xv, yv);
    ny = size(Zc,1); nx = size(Zc,2);

    % Trim axes
    xv = xv(:).'; yv = yv(:).';
    if numel(xv) ~= nx, xv = xv(1:nx); end
    if numel(yv) ~= ny, yv = yv(1:ny); end

    % Orientation / flips
    flipNotes = {};
    if xv(1) > xv(end)
        xv = fliplr(xv); Zc = fliplr(Zc);
        flipNotes{end+1} = 'flipped X';
    end
    if yv(1) > yv(end)
        yv = fliplr(yv); Zc = flipud(Zc);
        flipNotes{end+1} = 'flipped Y';
    end
    flipNoteStr = '';
    if ~isempty(flipNotes)
        flipNoteStr = [' (' strjoin(flipNotes, ', ') ')'];
    end

    % Downsample
    maxN = 1200;
    ix = round(linspace(1,nx,min(nx,maxN)));
    iy = round(linspace(1,ny,min(ny,maxN)));
    xv2 = xv(ix); yv2 = yv(iy); Z2 = Zc(iy,ix);

    % Project to EPSG:3031
    [xps,yps] = to3031(xv2, yv2, assume_deg, proj);

    % Plot
    figure('Color','w');
    imagesc([min(xps(:)) max(xps(:))],[min(yps(:)) max(yps(:))],Z2);
    set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); colorbar;
    title(sprintf('%s (SOURCE → EPSG:3031)%s\n[min=%.3g max=%.3g finite=%d/%d]', ...
        titleStr, flipNoteStr, ...
        min(Z(:),[],'omitnan'), max(Z(:),[],'omitnan'), ...
        nnz(isfinite(Z)), numel(Z)));
    xlabel('X (m)'); ylabel('Y (m)');
    try exportgraphics(gcf,out_png,'Resolution',180); catch, end
end

function qa_map_scatter_source(x, y, z, assume_deg, proj, titleStr, out_png)
    Mfin = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(Mfin); y = y(Mfin); z = z(Mfin);
    N = numel(x); maxN = 250000;
    if N > maxN, sel = randperm(N, maxN); x = x(sel); y = y(sel); z = z(sel); end
    if assume_deg, [xps, yps] = projfwd(proj, y, x); else, xps = x; yps = y; end
    figure('Color','w'); scatter(xps, yps, 6, z, 'filled'); axis equal tight; colorbar;
    title(sprintf('%s (SOURCE native → EPSG:3031)\n[min=%.3g max=%.3g finite=%d]', ...
        titleStr, min(z,[],'omitnan'), max(z,[],'omitnan'), numel(z)));
    xlabel('X (m)'); ylabel('Y (m)'); try exportgraphics(gcf, out_png, 'Resolution', 220); catch, end
end

%% ------------------ Stål (EPSG:3031; Q,U in W m^-2) ------------------
if doModels.Stal
    stal_file = 'data_ipics/aq1_01_20.nc';
    HF_stal   = ncread(stal_file,'Q');        % W/m^2
    UNC_stal  = ncread(stal_file,'U');        % W/m^2 (σ)
    x_stal    = ncread(stal_file,'X');        % meters (EPSG:3031)
    y_stal    = ncread(stal_file,'Y');
    
    HF_stal  = HF_stal * 1000;
    HF_stal = flipud(HF_stal);
    HF_stal = rot90(HF_stal, -1);
    
    UNC_stal = UNC_stal * 1000;
    UNC_stal = flipud(UNC_stal);
    UNC_stal = rot90(UNC_stal, -1);
    
    % Auto-detect axes (should be meters)
    [XS_stal, YS_stal, decided_stal] = to3031_auto(stal_file, 'X', 'Y', x_stal, y_stal, proj);
    qa_map_rect_source(x_stal, y_stal, HF_stal, strcmp(decided_stal,'degrees'), proj, ...
        'GHF Stal (SOURCE)', fullfile(qa_outdir,'SRC_GHF_Stal.png'));
    qa_map_rect_source(x_stal, y_stal, UNC_stal, strcmp(decided_stal,'degrees'), proj, ...
        'UNC Stal σ (SOURCE)', fullfile(qa_outdir,'SRC_UNC_Stal.png'));
    
    [F_stal, keptS, fracS, nfinS]   = make_interpolant(XS_stal(:), YS_stal(:), HF_stal(:));
    [U_stal, keptSU, fracSU, nfinSU]= make_interpolant(XS_stal(:), YS_stal(:), UNC_stal(:));
    fprintf('[Stal] mean kept=%d/%.0f%% finite=%d | unc kept=%d/%.0f%% finite=%d\n', ...
        keptS, 100*fracS, nfinS, keptSU, 100*fracSU, nfinSU);
    coverage_check('Stal', XS_stal, YS_stal, Xgrid, Ygrid);
    
    GHF_Stal_interp = F_stal(Xgrid, Ygrid);
    UNC_Stal_interp = U_stal(Xgrid, Ygrid);
    qa_map(GHF_Stal_interp, Xgrid, Ygrid, 'GHF Stal (mW m^{-2})', fullfile(qa_outdir,'map_GHF_Stal.png'));
    qa_hist(GHF_Stal_interp, 'GHF Stal', fullfile(qa_outdir,'hist_GHF_Stal.png'));
    qa_map(UNC_Stal_interp, Xgrid, Ygrid, 'UNC Stal (σ, mW m^{-2})', fullfile(qa_outdir,'map_UNC_Stal.png'));
    qa_hist(UNC_Stal_interp, 'UNC Stal (σ)', fullfile(qa_outdir,'hist_UNC_Stal.png'));
end

%% ------------------ Hazzard mean ------------------
if doModels.Hazzard
    haz_mean_file = 'data_ipics/HR24_GHF_mean_PS.grd';
    HF_haz = ncread(haz_mean_file,'z');
    x_haz = ncread(haz_mean_file,'x');  
    y_haz = ncread(haz_mean_file,'y');   
    
    HF_haz = to_mW_if_W(haz_mean_file,'z',HF_haz);
    HF_haz = clean_grid(HF_haz, haz_mean_file, 'z');
    HF_haz = flipud(HF_haz);
    HF_haz = rot90(HF_haz, -1);
    
    [x_hazg, y_hazg] = meshgrid(x_haz*1000, y_haz*1000);
    
    [F_haz, ~, ~] = make_interpolant(x_hazg(:), y_hazg(:), HF_haz(:));
    GHF_Hazzard_interp = F_haz(Xgrid, Ygrid);
    
    % === Std: already EPSG:3031 (km) ===
    haz_std_file = 'data_ipics/HR24_GHF_std_PS.grd';
    SD_haz = ncread(haz_std_file,'z');
    
    SD_haz = flipud(SD_haz);
    SD_haz = rot90(SD_haz, -1);

    x_haz2 = ncread(haz_std_file,'x');   % km
    y_haz2 = ncread(haz_std_file,'y');   % km
    
    % convert to meters for consistency
    [x_haz2g, y_haz2g] = meshgrid(x_haz2*1000, y_haz2*1000);
    
    SD_haz = to_mW_if_W(haz_std_file,'z',SD_haz);
    SD_haz = clean_grid(SD_haz, haz_std_file, 'z');
    
    [F_haz_std, kept_std, frac_std] = make_interpolant(x_haz2g(:), y_haz2g(:), SD_haz(:));
    UNC_Hazzard_interp = F_haz_std(Xgrid, Ygrid);
    
    fprintf('[Hazzard] mean=lat/lon projected, std=EPSG:3031 (km→m), kept=%d (%.2f%%)\n', ...
        kept_std, 100*frac_std);
    
    % QA plots
    qa_map(GHF_Hazzard_interp, Xgrid, Ygrid, 'GHF Hazzard (mW m^{-2})', fullfile(qa_outdir,'map_GHF_Hazzard.png'));
    qa_hist(GHF_Hazzard_interp, 'GHF Hazzard', fullfile(qa_outdir,'hist_GHF_Hazzard.png'));
    qa_map(UNC_Hazzard_interp, Xgrid, Ygrid, 'UNC Hazzard (σ, mW m^{-2})', fullfile(qa_outdir,'map_UNC_Hazzard.png'));
    qa_hist(UNC_Hazzard_interp, 'UNC Hazzard (σ)', fullfile(qa_outdir,'hist_UNC_Hazzard.png'));
end

%% ------------------ An ------------------
if doModels.An
    an_file = 'data_ipics/AN1-HF.grd';
    HF_an = ncread(an_file,'z');
    x_an  = ncread(an_file,'lon'); y_an = ncread(an_file,'lat');
    HF_an = to_mW_if_W(an_file,'z',HF_an);

    % Fix orientation if dimensions don't match
    if size(HF_an,1) == numel(x_an) && size(HF_an,2) == numel(y_an)
        HF_an = HF_an.';  % transpose
    end

    if y_an(1) > y_an(end)
        y_an = flipud(y_an);
        HF_an = flipud(HF_an);
    end

    [lon_an, lat_an] = meshgrid(x_an, y_an);
    
    % Project to EPSG:3031 (expects lat, lon as same-sized arrays)
    [Xps_an, Yps_an] = projfwd(proj, lat_an, lon_an);
    
    % Interpolant to your target grid
    [F_an, ~, ~] = make_interpolant(Xps_an(:), Yps_an(:), HF_an(:));
    GHF_An_interp = F_an(Xgrid, Ygrid);

    qa_map(GHF_An_interp, Xgrid, Ygrid, 'GHF An (mW m^{-2})', fullfile(qa_outdir,'map_GHF_An.png'));
    qa_hist(GHF_An_interp, 'GHF An', fullfile(qa_outdir,'hist_GHF_An.png'));
end   

%% ------------------ Fox Maule ------------------
if doModels.FoxMaule
    fox_file = 'data_ipics/fox_maule-hfmag.grd';
    HF_fox = ncread(fox_file,'z');
    x_fox  = ncread(fox_file,'x'); y_fox = ncread(fox_file,'y');
    HF_fox = to_mW_if_W(fox_file,'z',HF_fox);

    HF_fox = flipud(HF_fox);
    HF_fox = rot90(HF_fox, -1);
    
    [X_fox_ps, Y_fox_ps, decided_fox] = to3031_auto(fox_file, 'x','y', x_fox, y_fox, proj);
    qa_map_rect_source(x_fox, y_fox, HF_fox, strcmp(decided_fox,'degrees'), proj, ...
        'GHF FoxMaule (SOURCE)', fullfile(qa_outdir,'SRC_GHF_FoxMaule.png'));
    
    [F_fox, keptF, fracF, nfinF] = make_interpolant(X_fox_ps(:), Y_fox_ps(:), HF_fox(:));
    fprintf('[FoxMaule] kept=%d/%.0f%% finite=%d (axes=%s)\n', keptF, 100*fracF, nfinF, decided_fox);
    coverage_check('FoxMaule', X_fox_ps, Y_fox_ps, Xgrid, Ygrid);
    
    GHF_FoxMaule_interp = F_fox(Xgrid, Ygrid);
    qa_map(GHF_FoxMaule_interp, Xgrid, Ygrid, 'GHF FoxMaule (mW m^{-2})', fullfile(qa_outdir,'map_GHF_FoxMaule.png'));
    qa_hist(GHF_FoxMaule_interp, 'GHF FoxMaule', fullfile(qa_outdir,'hist_GHF_FoxMaule.png'));
end

%% ------------------ Martos (XYZ, EPSG:3031 assumed) ------------------
if doModels.Martos
    martos = readmatrix('data_ipics/Martos_GHF.xyz', 'FileType', 'text');
    xm = martos(:,1); ym = martos(:,2); zm = martos(:,3);
    qa_map_scatter_source(xm, ym, zm, false, proj, 'GHF Martos (SOURCE)', fullfile(qa_outdir,'SRC_GHF_Martos.png'));
    [F_martos, keptM, fracM, nfinM] = make_interpolant(xm, ym, zm);
    fprintf('[Martos] kept=%d/%.0f%% finite=%d\n', keptM, 100*fracM, nfinM);
    coverage_check('Martos mean', xm, ym, Xgrid, Ygrid);
    
    GHF_Martos_interp = F_martos(Xgrid, Ygrid);
    qa_map(GHF_Martos_interp, Xgrid, Ygrid, 'GHF Martos (mW m^{-2})', fullfile(qa_outdir,'map_GHF_Martos.png'));
    qa_hist(GHF_Martos_interp, 'GHF Martos', fullfile(qa_outdir,'hist_GHF_Martos.png'));
    
    martos_unc = readmatrix('data_ipics/Antarctic_GHF_uncertainty.xyz', 'FileType', 'text');
    xmu = martos_unc(:,1); ymu = martos_unc(:,2); zu = martos_unc(:,3);
    qa_map_scatter_source(xmu, ymu, zu, false, proj, 'UNC Martos σ (SOURCE)', fullfile(qa_outdir,'SRC_UNC_Martos.png'));
    [U_martos, keptMU, fracMU, nfinMU] = make_interpolant(xmu, ymu, zu);
    fprintf('[Martos σ] kept=%d/%.0f%% finite=%d\n', keptMU, 100*fracMU, nfinMU);
    coverage_check('Martos σ', xmu, ymu, Xgrid, Ygrid);
    
    UNC_Martos_interp = U_martos(Xgrid, Ygrid);
    qa_map(UNC_Martos_interp, Xgrid, Ygrid, 'UNC Martos (σ, mW m^{-2})', fullfile(qa_outdir,'map_UNC_Martos.png'));
    qa_hist(UNC_Martos_interp, 'UNC Martos (σ)', fullfile(qa_outdir,'hist_UNC_Martos.png'));
end

%% ------------------ Shen (XYZ lon/lat) ------------------
if doModels.Shen
    shen = readmatrix('data_ipics/shen.hf.v1.xyz', 'FileType', 'text');
    lonS = shen(:,1); latS = shen(:,2); zS = shen(:,3);
    qa_map_scatter_source(lonS, latS, zS, true, proj, 'GHF Shen (SOURCE)', fullfile(qa_outdir,'SRC_GHF_Shen.png'));
    [xpsS, ypsS] = projfwd(proj, latS, lonS);
    [F_shen, keptSh, fracSh, nfinSh] = make_interpolant(xpsS, ypsS, zS);
    fprintf('[Shen] kept=%d/%.0f%% finite=%d\n', keptSh, 100*fracSh, nfinSh);
    coverage_check('Shen mean', xpsS, ypsS, Xgrid, Ygrid);
    
    GHF_Shen_interp = F_shen(Xgrid, Ygrid);
    qa_map(GHF_Shen_interp, Xgrid, Ygrid, 'GHF Shen (mW m^{-2})', fullfile(qa_outdir,'map_GHF_Shen.png'));
    qa_hist(GHF_Shen_interp, 'GHF Shen', fullfile(qa_outdir,'hist_GHF_Shen.png'));
    
    shen_std = readmatrix('data_ipics/shen.unhf.v1.xyz', 'FileType','text');
    lonSu = shen_std(:,1); latSu = shen_std(:,2); zSu = shen_std(:,3);
    qa_map_scatter_source(lonSu, latSu, zSu, true, proj, 'UNC Shen σ (SOURCE)', fullfile(qa_outdir,'SRC_UNC_Shen.png'));
    [xpsSu, ypsSu] = projfwd(proj, latSu, lonSu);
    [STD_shen, keptShU, fracShU, nfinShU] = make_interpolant(xpsSu, ypsSu, zSu);
    fprintf('[Shen σ] kept=%d/%.0f%% finite=%d\n', keptShU, 100*fracShU, nfinShU);
    coverage_check('Shen σ', xpsSu, ypsSu, Xgrid, Ygrid);
    
    UNC_Shen_interp = STD_shen(Xgrid, Ygrid);
    qa_map(UNC_Shen_interp, Xgrid, Ygrid, 'UNC Shen (σ, mW m^{-2})', fullfile(qa_outdir,'map_UNC_Shen.png'));
    qa_hist(UNC_Shen_interp, 'UNC Shen (σ)', fullfile(qa_outdir,'hist_UNC_Shen.png'));
end

%% ------------------ Losing (CSV lon/lat mean + min/max) ------------------
if doModels.Losing
    losing_tbl = readtable('data_ipics/HF_Min_Max_MaxAbs-1.csv');
    lonL = losing_tbl{:,1}; latL = losing_tbl{:,2};
    zL   = losing_tbl{:,3}; zLmin = losing_tbl{:,4}; zLmax = losing_tbl{:,5};
    
    qa_map_scatter_source(lonL, latL, zL, true, proj, 'GHF Losing (SOURCE)',  fullfile(qa_outdir,'SRC_GHF_Losing.png'));
    qa_map_scatter_source(lonL, latL, zLmin, true, proj, 'BMIN Losing (SOURCE)', fullfile(qa_outdir,'SRC_BMIN_Losing.png'));
    qa_map_scatter_source(lonL, latL, zLmax, true, proj, 'BMAX Losing (SOURCE)', fullfile(qa_outdir,'SRC_BMAX_Losing.png'));
    
    [xpsL, ypsL] = projfwd(proj, latL, lonL);
    [F_L, keptL, fracL, nfinL] = make_interpolant(xpsL, ypsL, zL);
    [F_Lmin, keptLmin, fracLmin, nfinLmin] = make_interpolant(xpsL, ypsL, zLmin);
    [F_Lmax, keptLmax, fracLmax, nfinLmax] = make_interpolant(xpsL, ypsL, zLmax);
    fprintf('[Losing] mean kept=%d/%.0f%%; min kept=%d/%.0f%%; max kept=%d/%.0f%%\n', ...
        keptL, 100*fracL, keptLmin, 100*fracLmin, keptLmax, 100*fracLmax);
    coverage_check('Losing mean', xpsL, ypsL, Xgrid, Ygrid);
    
    GHF_Losing_interp  = F_L(Xgrid, Ygrid);
    BMIN_Losing_interp = F_Lmin(Xgrid, Ygrid);
    BMAX_Losing_interp = F_Lmax(Xgrid, Ygrid);
    qa_map(GHF_Losing_interp,  Xgrid, Ygrid, 'GHF Losing (mW m^{-2})', fullfile(qa_outdir,'map_GHF_Losing.png'));
    qa_hist(GHF_Losing_interp, 'GHF Losing', fullfile(qa_outdir,'hist_GHF_Losing.png'));
    qa_map(BMIN_Losing_interp, Xgrid, Ygrid, 'Losing MIN (mW m^{-2})', fullfile(qa_outdir,'map_BMIN_Losing.png'));
    qa_map(BMAX_Losing_interp, Xgrid, Ygrid, 'Losing MAX (mW m^{-2})', fullfile(qa_outdir,'map_BMAX_Losing.png'));
end

%% ------------------ Save everything to datasets_for_gmin/ (+ GeoTIFFs) ------------------
if doModels.save
    outdir = 'datasets_for_gmin/';
    if ~exist(outdir,'dir'), mkdir(outdir); end

    %% Means (MAT + GeoTIFF)
    if exist('GHF_Stal_interp','var')
        save(fullfile(outdir,'GHF_Stal_interp.mat'), 'GHF_Stal_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_Stal_interp.tif'), ...
            GHF_Stal_interp, Xgrid, Ygrid);
    end

    if exist('GHF_Hazzard_interp','var')
        save(fullfile(outdir,'GHF_Hazzard_interp.mat'), 'GHF_Hazzard_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_Hazzard_interp.tif'), ...
            GHF_Hazzard_interp, Xgrid, Ygrid);
    end

    if exist('GHF_An_interp','var')
        save(fullfile(outdir,'GHF_An_interp.mat'), 'GHF_An_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_An_interp.tif'), ...
            GHF_An_interp, Xgrid, Ygrid);
    end

    if exist('GHF_Martos_interp','var')
        save(fullfile(outdir,'GHF_Martos_interp.mat'), 'GHF_Martos_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_Martos_interp.tif'), ...
            GHF_Martos_interp, Xgrid, Ygrid);
    end

    if exist('GHF_Shen_interp','var')
        save(fullfile(outdir,'GHF_Shen_interp.mat'), 'GHF_Shen_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_Shen_interp.tif'), ...
            GHF_Shen_interp, Xgrid, Ygrid);
    end

    if exist('GHF_Losing_interp','var')
        save(fullfile(outdir,'GHF_Losing_interp.mat'), 'GHF_Losing_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_Losing_interp.tif'), ...
            GHF_Losing_interp, Xgrid, Ygrid);
    end

    if exist('GHF_FoxMaule_interp','var')
        save(fullfile(outdir,'GHF_FoxMaule_interp.mat'), 'GHF_FoxMaule_interp');
        write_geotiff_epsg3031(fullfile(outdir,'GHF_FoxMaule_interp.tif'), ...
            GHF_FoxMaule_interp, Xgrid, Ygrid);
    end

    %% Uncertainties / bounds (MAT + GeoTIFF)
    if exist('UNC_Stal_interp','var')
        save(fullfile(outdir,'UNC_Stal_interp.mat'), 'UNC_Stal_interp');
        write_geotiff_epsg3031(fullfile(outdir,'UNC_Stal_interp.tif'), ...
            UNC_Stal_interp, Xgrid, Ygrid);
    end

    if exist('UNC_Hazzard_interp','var')
        save(fullfile(outdir,'UNC_Hazzard_interp.mat'), 'UNC_Hazzard_interp');
        write_geotiff_epsg3031(fullfile(outdir,'UNC_Hazzard_interp.tif'), ...
            UNC_Hazzard_interp, Xgrid, Ygrid);
    end

    if exist('UNC_Martos_interp','var')
        save(fullfile(outdir,'UNC_Martos_interp.mat'), 'UNC_Martos_interp');
        write_geotiff_epsg3031(fullfile(outdir,'UNC_Martos_interp.tif'), ...
            UNC_Martos_interp, Xgrid, Ygrid);
    end

    if exist('UNC_Shen_interp','var')
        save(fullfile(outdir,'UNC_Shen_interp.mat'), 'UNC_Shen_interp');
        write_geotiff_epsg3031(fullfile(outdir,'UNC_Shen_interp.tif'), ...
            UNC_Shen_interp, Xgrid, Ygrid);
    end

    if exist('BMIN_Losing_interp','var')
        save(fullfile(outdir,'BMIN_Losing_interp.mat'), 'BMIN_Losing_interp');
        write_geotiff_epsg3031(fullfile(outdir,'BMIN_Losing_interp.tif'), ...
            BMIN_Losing_interp, Xgrid, Ygrid);
    end

    if exist('BMAX_Losing_interp','var')
        save(fullfile(outdir,'BMAX_Losing_interp.mat'), 'BMAX_Losing_interp');
        write_geotiff_epsg3031(fullfile(outdir,'BMAX_Losing_interp.tif'), ...
            BMAX_Losing_interp, Xgrid, Ygrid);
    end

    fprintf('All interpolated grids saved as MAT + GeoTIFF in %s\n', outdir);
    fprintf('All QA maps (source + target) saved in %s\n', qa_outdir);
end

function write_geotiff_epsg3031(filename, Z, Xgrid, Ygrid, nodataVal)
% WRITE_GEOTIFF_EPSG3031  Save a grid on Xgrid/Ygrid as a GeoTIFF (EPSG:3031).
%
%   - Assumes Xgrid, Ygrid are cell centers in meters (EPSG:3031).
%   - Expands world limits by half a cell to get cell edges.
%   - Replaces non-finite values with nodataVal (default -9999).

    if nargin < 5 || isempty(nodataVal)
        nodataVal = -9999;
    end

    % Replace NaN/Inf with nodata
    Zout = Z;
    bad = ~isfinite(Zout);
    Zout(bad) = nodataVal;

    % Assume regular grid; use first row/column as 1D axes
    xv = Xgrid(1, :);
    yv = Ygrid(:, 1);

    dx = median(diff(xv), 'omitnan');
    dy = median(diff(yv), 'omitnan');

    % Cell-edge world limits (half-cell padding)
    xmin = min(xv) - dx/2;
    xmax = max(xv) + dx/2;
    ymin = min(yv) - dy/2;
    ymax = max(yv) + dy/2;

    % Build spatial referencing object
    R = maprefcells([xmin xmax], [ymin ymax], size(Zout));

    % Try to make rows start from south (y increasing upward),
    % but don't die if this MATLAB version doesn't like it.
    try
        R.RowsStartFrom = 'south';
    catch
        % leave default
    end

    % DO NOT touch ColumnsStartFrom here — some versions are picky.
    % Just use whatever default maprefcells gave us.

    % Write GeoTIFF (EPSG:3031)
    geotiffwrite(filename, Zout, R, 'CoordRefSysCode', 3031);
    fprintf('[GeoTIFF] wrote %s (EPSG:3031, nodata=%g)\n', filename, nodataVal);
end

