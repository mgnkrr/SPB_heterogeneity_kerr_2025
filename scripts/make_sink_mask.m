function make_sink_mask()
% MAKE_SINK_MASK
% Build sink mask from hydraulic head on a 2 m grid (REMA-based), then
% regrid the resulting sink mask to the 100 m analysis grid used by the
% GHF evaluation / bootstrap code.
%
% Outputs (on 100 m grid):
%   datasets_for_gmin/sink_mask_new.mat  -> sink_mask, Xgrid, Ygrid, cfg
%   datasets_for_gmin/sink_mask_comp.mat -> sink_mask_comp, comp_id, Xgrid, Ygrid
%
% Optional outputs (specularity vs sinks):
%   datasets_for_gmin/spec_sink_stats/spec_sink_stats.mat
%   datasets_for_gmin/spec_sink_stats/png/*.png
%   datasets_for_gmin/spec_sink_stats/wet_masks/*.mat   (optional; compact index form)
%
% Notes:
% - Sink analysis is done on *sink pixels* (sink_mask), not components.
% - Specularity convention: Q==9999 treated as missing (NaN).

%% -------------------- CONFIG --------------------
conn        = 8;                      % 4 or 8 connectivity (for 100 m comps)
min_size    = 5;                      % min pixels per component
out_sinks   = fullfile('datasets_for_gmin','sink_mask_new.mat');
out_comp    = fullfile('datasets_for_gmin','sink_mask_comp.mat');

addpath('/Users/megankerr/Documents/Location-oldest-ice/data_ipics/Bed Head v1.0/Bed Head');

%% -------------------- Specularity vs sinks analysis --------------------
doSpecAnalysis   = true;

spec_file        = fullfile('datasets_for_gmin','specularity.mat');
spec_outdir      = fullfile('datasets_for_gmin','spec_sink_stats');

wet_thr0         = 0.2;     % reference threshold to report/mark on plots
nSweep           = 80;        % number of sweep thresholds (quantile-based)
minFiniteFrac    = 0.02;      % keep thresholds where wetFrac in [minFiniteFrac, 1-minFiniteFrac]

% Save per-threshold wet masks?
% IMPORTANT: saving full logical grids is huge; we save compact linear indices instead.
saveSweepMasks   = true;
maskFormat       = "index";   % "index" (recommended) or "full" (HUGE)

%% -------------------- Load 2 m rasters --------------------
disp('Loading 2 m input data...');
load('datasets_for_gmin/coldex_bedelv_2m.mat',  'B_2m', 'X2m', 'Y2m');
load('datasets_for_gmin/REMA_srfelv_2m.mat',    'S_rem_2m');
load('datasets_for_gmin/coldex_srfelv_2m.mat',  'S_coldex_2m');

assert(isequal(size(B_2m), size(S_rem_2m), size(S_coldex_2m), size(X2m), size(Y2m)), ...
    'Size mismatch among 2 m rasters');

% Basic diagnostics on 2 m grid
xv2 = X2m(1,:);
yv2 = Y2m(:,1);
dx2 = median(diff(xv2));
dy2 = median(diff(yv2));
fprintf('[2m grid] size=%dx%d | X %s (dx=%.3f m) | Y %s (dy=%.3f m)\n', ...
    size(B_2m,1), size(B_2m,2), tern(all(diff(xv2)>0),'asc','desc'), dx2, tern(all(diff(yv2)>0),'asc','desc'), dy2);

%% -------------------- Ice thickness on 2 m grid --------------------
proj = projcrs(3031);
[lat2, ~] = projinv(proj, X2m, Y2m);

mask_coldex = (lat2 <= -88);    % use COLDEX surface near the pole, REMA elsewhere
combined_srfelv_2m = S_rem_2m;
combined_srfelv_2m(mask_coldex) = S_coldex_2m(mask_coldex);

% Area where COLDEX replaces REMA
valid_coldex = mask_coldex & isfinite(S_coldex_2m);
pixel_area_m2   = abs(dx2 * dy2);
n_coldex_pixels = nnz(valid_coldex);
area_coldex_km2 = n_coldex_pixels * pixel_area_m2 / 1e6;
n_finite = nnz(isfinite(combined_srfelv_2m));
frac_coldex_pct = 100 * n_coldex_pixels / max(n_finite,1);

fprintf('[COLDEX replace] %d px -> %.2f km^2 (%.2f%% of finite 2 m surface grid)\n', ...
    n_coldex_pixels, area_coldex_km2, frac_coldex_pct);

cfg = struct();
cfg.coldex_replace_npix  = n_coldex_pixels;
cfg.coldex_replace_km2   = area_coldex_km2;
cfg.coldex_replace_pct   = frac_coldex_pct;

% Ice thickness (nonnegative)
thk_2m = max(combined_srfelv_2m - B_2m, 0);

%% -------------------- Normalize to north-up for GRIDobj (2 m) --------
% GRIDobj expects y strictly decreasing down rows.
flipY2 = (yv2(1) < yv2(end));  % true if south-up (ascending y)
if flipY2
    disp('[2m orient] Detected south-up input (Y ascending). Flipping to north-up for processing.');
    Bz2       = flipud(B_2m);
    Thkz2     = flipud(thk_2m);
    yv2_north = flipud(yv2);
else
    Bz2       = B_2m;
    Thkz2     = thk_2m;
    yv2_north = yv2;
end
xv2_east = xv2;  % assume X increases eastward in EPSG:3031

%% -------------------- GRIDobjs for TopoToolbox (2 m) -----------------
disp('Creating 2 m GRID objects (north-up) ...');
bed_dem_2m = GRIDobj(xv2_east, yv2_north, Bz2);
thk_dem_2m = GRIDobj(xv2_east, yv2_north, Thkz2);

%% -------------------- Hydraulic head on 2 m grid ---------------------
try
    disp('Calculating hydraulic head with bedhead() on 2 m grid ...');
    head_dem_2m = bedhead('bed', bed_dem_2m, 'thickness', thk_dem_2m);
catch
    warning('bedhead() not found; using bed elevation as head proxy (2 m).');
    head_dem_2m = bed_dem_2m;
end

%% -------------------- Fill sinks on 2 m grid -------------------------
disp('Computing DEM fill (2 m, north-up) ...');
DEMf_2m = fillsinks(head_dem_2m);

sinks_depth_north_2m = DEMf_2m.Z - head_dem_2m.Z;
sink_mask_north_2m   = (sinks_depth_north_2m > 0) & isfinite(sinks_depth_north_2m);

% Flip back to original orientation (2 m grid)
if flipY2
    sink_mask_2m = flipud(sink_mask_north_2m);
else
    sink_mask_2m = sink_mask_north_2m;
end

%% -------------------- Regrid sinks to 100 m analysis grid ------------
disp('Regridding 2 m sink mask to 100 m analysis grid ...');
L100 = load('datasets_for_gmin/coldex_bedelv.mat', 'B', 'Xgrid', 'Ygrid');
Xgrid = L100.Xgrid;
Ygrid = L100.Ygrid;

% Ensure ascending axes for interpolant inputs
xv2a = xv2; yv2a = yv2; sink2a = sink_mask_2m;
if numel(xv2a)>1 && xv2a(2) < xv2a(1)
    xv2a = fliplr(xv2a);
    sink2a = fliplr(sink2a);
end
if numel(yv2a)>1 && yv2a(2) < yv2a(1)
    yv2a = flipud(yv2a);
    sink2a = flipud(sink2a);
end

% Sample at 100 m grid centers
F_sink = griddedInterpolant({yv2a, xv2a}, double(sink2a), 'nearest', 'none');
sink_interp = F_sink(Ygrid, Xgrid);
sink_mask = isfinite(sink_interp) & (sink_interp > 0.5);

fprintf('[100m sinks] grid size=%dx%d | sink px=%d (%.3f%%)\n', ...
    size(sink_mask,1), size(sink_mask,2), nnz(sink_mask), 100*nnz(sink_mask)/numel(sink_mask));

%% -------------------- Save sink pixel mask (100 m) --------------------
if ~exist(fileparts(out_sinks), 'dir'), mkdir(fileparts(out_sinks)); end
save(out_sinks, 'sink_mask', 'Xgrid', 'Ygrid', 'cfg', '-v7.3');
fprintf('Saved 100 m sink mask to: %s\n', out_sinks);

%% -------------------- Connected components on 100 m grid -------------
disp('Finding connected components on 100 m grid ...');
CC = bwconncomp(sink_mask, conn);
fprintf('  initial components: %d\n', CC.NumObjects);

keepC = cellfun(@numel, CC.PixelIdxList) >= min_size;
sink_mask_comp = false(size(sink_mask));
if any(keepC)
    sink_mask_comp(vertcat(CC.PixelIdxList{keepC})) = true;
end

CC2 = bwconncomp(sink_mask_comp, conn);
comp_id = uint32(labelmatrix(CC2));
fprintf('  kept components: %d (removed %d)\n', CC2.NumObjects, CC.NumObjects - CC2.NumObjects);

%% -------------------- Save component labels (100 m) -------------------
if ~exist(fileparts(out_comp), 'dir'), mkdir(fileparts(out_comp)); end
save(out_comp, 'sink_mask_comp', 'comp_id', 'Xgrid', 'Ygrid', '-v7');
fprintf('Saved 100 m component mask/labels -> %s\n', out_comp);

%% -------------------- Specularity vs sinks (hist + sweep) ------------
if doSpecAnalysis
    try
        analyze_specularity_vs_sinks_pixels( ...
            spec_file, spec_outdir, sink_mask, Xgrid, Ygrid, wet_thr0, nSweep, minFiniteFrac, ...
            saveSweepMasks, maskFormat);
    catch ME
        fprintf(2,'Specularity analysis failed: %s\n', ME.message);
        % disp(getReport(ME,'extended'));
    end
end

%% -------------------- Quicklook (sink pixels over 100 m bed) ---------
figure('Color','w','Name','Sink pixels over 100 m bed');
imagesc(Xgrid(1,:), Ygrid(:,1), L100.B);
axis image; set(gca,'YDir','normal'); colormap gray; hold on
title(sprintf('Sink pixels (100 m) | sink px=%d (%.2f%%)', nnz(sink_mask), 100*nnz(sink_mask)/numel(sink_mask)));
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on
contour(Xgrid, Ygrid, double(sink_mask), [0.5 0.5], 'm', 'LineWidth', 1.5);

end

%% ======================= Local helper functions =======================

function s = tern(cond, a, b)
if all(cond), s = a; else, s = b; end
end

function analyze_specularity_vs_sinks_pixels(spec_file, outdir, sink_mask, Xgrid, Ygrid, wet_thr0, nSweep, minFiniteFrac, saveMasks, maskFormat)
% Robust analysis of specularity Q vs sink pixels (sink_mask).
%
% Outputs:
%   outdir/spec_sink_stats.mat
%   outdir/png/*.png
%   outdir/wet_masks/*.mat (optional)

    if ~exist(outdir,'dir'), mkdir(outdir); end
    pngdir  = fullfile(outdir,'png');       if ~exist(pngdir,'dir'), mkdir(pngdir); end
    maskdir = fullfile(outdir,'wet_masks'); if saveMasks && ~exist(maskdir,'dir'), mkdir(maskdir); end

    assert(isfile(spec_file), 'Specularity file not found: %s', spec_file);
    L = load(spec_file, 'Q', 'Xgrid', 'Ygrid');
    Q = double(L.Q);

    % --- Grid sanity ---
    assert(isequal(size(Q), size(sink_mask)), ...
        'Q size %s != sink_mask size %s', mat2str(size(Q)), mat2str(size(sink_mask)));
    assert(isequal(size(L.Xgrid), size(Xgrid)) && isequal(size(L.Ygrid), size(Ygrid)), ...
        'Specularity Xgrid/Ygrid do not match sink grid. Regrid specularity first.');

    % --- Clean missing / nonphysical ---
    %Q(Q==9999) = NaN;
    Q(Q<=0) = NaN;  % comment out if zeros are meaningful

    % --- Valid domain ---
    valid = isfinite(Q);
    nV = nnz(valid);
    fprintf('[spec-sink] valid pixels: %d / %d (%.1f%%)\n', nV, numel(Q), 100*nV/numel(Q));
    if nV == 0, error('No finite Q values after cleaning.'); end

    sink    = logical(sink_mask) & valid;
    nonsink = (~logical(sink_mask)) & valid;

    qS  = Q(sink);
    qNS = Q(nonsink);

    n1 = numel(qS); n2 = numel(qNS);
    fprintf('[spec-sink] nSink(valid)=%d | nNonSink(valid)=%d\n', n1, n2);
    if n1 == 0 || n2 == 0, error('Need nonempty sink & nonsink Q samples.'); end

    %% -------------------- Tests (MW, KS) --------------------
    % MW ranksum
    try
        p_mw = ranksum(qS, qNS);
    catch
        nsub = 2e6;
        iS  = randsample(n1,  min(nsub, n1));
        iNS = randsample(n2,  min(nsub, n2));
        p_mw = ranksum(qS(iS), qNS(iNS));
    end

    % % KS
    try
        [~, p_ks, D] = kstest2(qS, qNS);
    catch
        nsub = 2e6;
        iS  = randsample(n1,  min(nsub, n1));
        iNS = randsample(n2,  min(nsub, n2));
        [~, p_ks, D] = kstest2(qS(iS), qNS(iNS));
    end

    %% -------------------- PS (paired Monte Carlo) + bootstrap CI --------------------
    ps_cfg = struct();
    ps_cfg.nPairs     = 2e6;   % point estimate
    ps_cfg.B          = 300;   % bootstrap replicates
    ps_cfg.nPairsBoot = 3e5;   % pairs per bootstrap
    ps_cfg.seed       = 1;

    [PS, PS_CI] = prob_superiority_mc(qS, qNS, ps_cfg);

    % Distribution summaries
    med_in  = median(qS,'omitnan');
    med_out = median(qNS,'omitnan');
    dmed    = med_in - med_out;
    frac_in  = mean(qS >= wet_thr0);
    frac_out = mean(qNS >= wet_thr0);

    fprintf('\nSpecularity vs sinks: statistical summary\n');
    fprintf('n(valid) = %d | n(in sinks) = %d | n(out sinks) = %d\n', nV, n1, n2);
    fprintf('Mann-Whitney ranksum p = %.3e (log10 p ≈ %.2f)\n', p_mw, log10(max(p_mw, realmin)));
    fprintf('KS test p = %.3e | D = %.3f\n', p_ks, D);
    fprintf('Probability of superiority PS = %.3f (bootstrap 95%% CI: [%.3f, %.3f])\n', PS, PS_CI(1), PS_CI(2));
    fprintf('Median(Q): in = %.3f | out = %.3f | Δ = %.3f\n', med_in, med_out, dmed);
    fprintf('Frac (@ >= %.3f): in = %.3f | out = %.3f\n', wet_thr0, frac_in, frac_out);

    if ~(PS_CI(1) <= PS && PS <= PS_CI(2))
        error('PS %.3f is not inside bootstrap CI [%.3f, %.3f]. Pairing/stale-var issue.', ...
            PS, PS_CI(1), PS_CI(2));
    end

    %% -------------------- Distribution plots (PDF / logPDF / CDF) --------------------
    pooled = Q(valid);
    lo = prctile(pooled, 0.5);
    hi = prctile(pooled, 99.5);
    nbin = 60;
    edges = linspace(lo, hi, nbin+1);
    qcent = 0.5*(edges(1:end-1) + edges(2:end));

    hs = histcounts(qS,  edges, 'Normalization','pdf');
    hn = histcounts(qNS, edges, 'Normalization','pdf');

    fig = figure('Color','w','Visible','off'); hold on; box on
    plot(qcent, hn, 'k-', 'LineWidth', 2);
    plot(qcent, hs, 'b-', 'LineWidth', 2);
    grid on
    xlabel('Specularity content (Q)');
    ylabel('Probability density');
    title('Specularity inside vs outside hydrostatic sinks (pixels)');
    legend({sprintf('outside sinks (n=%d)', n2), sprintf('inside sinks (n=%d)', n1)}, 'Location','best');
    xline(wet_thr0, '--', sprintf('Q = %.3f', wet_thr0), 'LineWidth', 1.2, 'Interpreter','none');
    exportgraphics(fig, fullfile(pngdir,'histpdf_Q_inSink_vs_outSink_PIXELS.png'), 'Resolution', 220);
    close(fig);

    fig = figure('Color','w','Visible','off'); hold on; box on
    plot(qcent, hn, 'k-', 'LineWidth', 2);
    plot(qcent, hs, 'b-', 'LineWidth', 2);
    set(gca,'YScale','log');
    grid on
    xlabel('Specularity content (Q)');
    ylabel('Probability density (log scale)');
    title('Specularity PDF (log-y) inside vs outside sinks (pixels)');
    legend({'outside sinks','inside sinks'}, 'Location','best');
    xline(wet_thr0, '--', sprintf('Q = %.3f', wet_thr0), 'LineWidth', 1.2, 'Interpreter','none');
    exportgraphics(fig, fullfile(pngdir,'histpdf_log_Q_inSink_vs_outSink_PIXELS.png'), 'Resolution', 220);
    close(fig);

    fig = figure('Color','w','Visible','off'); hold on; box on
    [F1,X1] = ecdf(qS);
    [F2,X2] = ecdf(qNS);
    plot(X2, F2, 'k-', 'LineWidth', 2);
    plot(X1, F1, 'b-', 'LineWidth', 2);
    grid on
    xlabel('Specularity content (Q)');
    ylabel('CDF');
    title('CDF of specularity inside vs outside sinks (pixels)');
    legend({sprintf('outside sinks (n=%d)', n2), sprintf('inside sinks (n=%d)', n1)}, 'Location','best');
    xline(wet_thr0, '--', sprintf('Q = %.3f', wet_thr0), 'LineWidth', 1.2, 'Interpreter','none');
    exportgraphics(fig, fullfile(pngdir,'cdf_Q_inSink_vs_outSink_PIXELS.png'), 'Resolution', 220);
    close(fig);

    %% -------------------- Quantile shift summary --------------------
    qLevels = [10 25 50 75 90 95 99];
    q_in  = prctile(qS,  qLevels);
    q_out = prctile(qNS, qLevels);
    dq    = q_in - q_out;

    %% -------------------- Robustness: stride subsampling --------------------
    robustness = struct();
    targetFrac = 0.02;
    stride_k = max(1, round(sqrt(1/targetFrac)));  % ~7 => ~2%
    stride_k = min(stride_k, 25);

    [Qsub, sinksub] = subsample_stride(Q, sink_mask, valid, stride_k);
    qS_sub  = Qsub(sinksub);
    qNS_sub = Qsub(~sinksub);

    if numel(qS_sub) >= 1000 && numel(qNS_sub) >= 1000
        p_mw_sub = ranksum(qS_sub, qNS_sub);
        [~, p_ks_sub, D_sub] = kstest2(qS_sub, qNS_sub);
        [~,~,mw_sub] = ranksum(qS_sub, qNS_sub);

        n1s = numel(qS_sub);
        n2s = numel(qNS_sub);
        U_sub = mw_sub.ranksum - n1s*(n1s+1)/2;
        PS_sub = U_sub / (n1s*n2s);

        %PS_rank_sub = mw_sub.ranksum / (numel(qS_sub)*numel(qNS_sub));

        robustness.did_run = true;
        robustness.stride_k = stride_k;
        robustness.nSink = numel(qS_sub);
        robustness.nNonSink = numel(qNS_sub);
        robustness.p_mw = p_mw_sub;
        robustness.p_ks = p_ks_sub;
        robustness.D = D_sub;
        robustness.PS_rank = PS_sub;

        fprintf('[spec-sink][robust] stride=%d => ranksum p=%.3g | KS p=%.3g (D=%.3f) | PS=%.3f\n', ...
            stride_k, p_mw_sub, p_ks_sub, D_sub, PS_sub);

        robustness.PS = PS_sub;

    else
        robustness.did_run = false;
        robustness.stride_k = stride_k;
        robustness.note = 'Not enough subsampled points for robust tests.';
        fprintf('[spec-sink][robust] skipped: too few subsampled points.\n');
    end

    %% -------------------- Threshold sweep (quantile-based) --------------------
    qs = pooled(:);
    pp = linspace(0.01, 0.99, nSweep);
    thr = quantile(qs, pp);
    thr = unique(thr(:), 'stable');

    qv = Q(valid);
    sinkV_full = sink(valid);  % sink on valid-domain vector

    wetFrac = nan(size(thr));
    for i = 1:numel(thr)
        wetFrac(i) = mean(qv >= thr(i));
    end

    keep = wetFrac >= minFiniteFrac & wetFrac <= (1-minFiniteFrac);
    thr = thr(keep);
    wetFrac = wetFrac(keep);

    S = struct();
    S.thr = thr;
    S.wetFrac = wetFrac;
    S.TPR = nan(size(thr));
    S.FPR = nan(size(thr));
    S.precision = nan(size(thr));
    S.specificity = nan(size(thr));
    S.balAcc = nan(size(thr));
    S.MCC = nan(size(thr));
    S.F1 = nan(size(thr));
    S.oddsRatio = nan(size(thr));

    if saveMasks && maskFormat=="index"
        full_idx_valid = find(valid);
    end

    for i = 1:numel(thr)
        wetV = (qv >= thr(i));
        Mi = confusion_metrics_vector(wetV, sinkV_full);

        S.TPR(i)         = Mi.TPR;
        S.FPR(i)         = Mi.FPR;
        S.precision(i)   = Mi.precision;
        S.specificity(i) = Mi.specificity;
        S.balAcc(i)      = Mi.balAcc;
        S.MCC(i)         = Mi.MCC;
        S.F1(i)          = Mi.F1;
        S.oddsRatio(i)   = Mi.oddsRatio;

        if saveMasks
            thr_i = thr(i);
            wetFrac_i = wetFrac(i);
            outname = fullfile(maskdir, sprintf('wetmask_thr_%03d.mat', i));

            if maskFormat=="index"
                wet_full_idx = uint32(full_idx_valid(wetV));
                save(outname, 'thr_i', 'wetFrac_i', 'wet_full_idx', '-v7.3');
            elseif maskFormat=="full"
                wet_full = false(size(Q));
                wet_full(valid) = wetV;
                save(outname, 'thr_i', 'wetFrac_i', 'wet_full', '-v7.3');
            else
                error('Unknown maskFormat: %s', maskFormat);
            end
        end
    end

    [~, k0] = min(abs(S.thr - wet_thr0));
    fprintf('[spec-sink] thr≈%.3f: MCC=%.3f | balAcc=%.3f | TPR=%.3f | FPR=%.3f | precision=%.3f\n', ...
        S.thr(k0), S.MCC(k0), S.balAcc(k0), S.TPR(k0), S.FPR(k0), S.precision(k0));

    [~, imcc] = max(S.MCC);
    [~, ij]   = max(S.TPR - S.FPR);
    best = struct();
    best.byMCC = struct('thr', S.thr(imcc), 'MCC', S.MCC(imcc), 'wetFrac', S.wetFrac(imcc));
    best.byJ   = struct('thr', S.thr(ij),   'J', (S.TPR(ij)-S.FPR(ij)), 'wetFrac', S.wetFrac(ij));

    fig = figure('Color','w','Visible','off'); box on
    plot(S.thr, S.MCC, 'LineWidth', 2); grid on
    xlabel('Threshold t (wet if Q >= t)', 'Interpreter','none');
    ylabel('MCC');
    title('Specularity threshold sweep vs sinks (pixels)');
    xline(wet_thr0, '--', sprintf('Q=%.3f', wet_thr0), 'LineWidth', 1.2, 'Interpreter','none');
    exportgraphics(fig, fullfile(pngdir,'sweep_MCC_vs_threshold_PIXELS.png'), 'Resolution', 220);
    close(fig);

    fig = figure('Color','w','Visible','off'); hold on; box on
    plot(S.thr, S.TPR, 'LineWidth', 2);
    plot(S.thr, S.FPR, 'LineWidth', 2);
    grid on
    xlabel('Threshold t (wet if Q >= t)', 'Interpreter','none');
    ylabel('Rate');
    title('TPR/FPR vs specularity threshold (pixels)');
    legend({'TPR (recall)','FPR'}, 'Location','best');
    xline(wet_thr0, '--', sprintf('Q=%.3f', wet_thr0), 'LineWidth', 1.2, 'Interpreter','none');
    exportgraphics(fig, fullfile(pngdir,'sweep_TPR_FPR_vs_threshold_PIXELS.png'), 'Resolution', 220);
    close(fig);

    fig = figure('Color','w','Visible','off'); hold on; box on
    plot(S.FPR, S.TPR, 'LineWidth', 2);
    plot([0 1],[0 1],'--');
    grid on
    xlabel('FPR'); ylabel('TPR');
    title('ROC implied by specularity threshold sweep (pixels)');
    exportgraphics(fig, fullfile(pngdir,'sweep_ROC_PIXELS.png'), 'Resolution', 220);
    close(fig);

    %% -------------------- Save stats bundle --------------------
    stats = struct();
    stats.wet_thr0 = wet_thr0;

    stats.nValid = nV;
    stats.nSinkValid = nnz(sink);
    stats.nNonSinkValid = nnz(nonsink);

    stats.tests = struct();
    stats.tests.mannwhitney = struct('p', p_mw, 'n_in', n1, 'n_out', n2);
    stats.tests.ks = struct('p', p_ks, 'D', D);

    stats.effect = struct();
    stats.effect.PS = PS;
    stats.effect.PS_CI = PS_CI;
    stats.effect.PS_cfg = ps_cfg;

    stats.dist = struct();
    stats.dist.median_in  = med_in;
    stats.dist.median_out = med_out;
    stats.dist.delta_median = dmed;
    stats.dist.fracWet_in  = frac_in;
    stats.dist.fracWet_out = frac_out;

    stats.quantiles = struct();
    stats.quantiles.levels = qLevels;
    stats.quantiles.in  = q_in;
    stats.quantiles.out = q_out;
    stats.quantiles.delta_in_minus_out = dq;

    stats.robustness = robustness;

    stats.sweep = S;
    stats.best = best;
    stats.saveSweepMasks = saveMasks;
    stats.maskFormat = maskFormat;

    save(fullfile(outdir,'spec_sink_stats.mat'), 'stats', '-v7.3');
    fprintf('[spec-sink] wrote: %s\n', fullfile(outdir,'spec_sink_stats.mat'));
    fprintf('[spec-sink] plots: %s\n', pngdir);
    if saveMasks
        fprintf('[spec-sink] wet masks: %s (format=%s)\n', maskdir, maskFormat);
    end

    %% -------------------- Extra summary figure (stats in one place) --------------------
    fig = figure('Color','w','Visible','off');
    axis off

    lines = {};
    lines{end+1} = sprintf('Specularity vs sinks: statistical summary');
    lines{end+1} = sprintf('n(valid) = %d | n(in sinks) = %d | n(out sinks) = %d', nV, n1, n2);
    lines{end+1} = sprintf('Mann–Whitney ranksum p = %.3g', p_mw);
    lines{end+1} = sprintf('KS test p = %.3g | D = %.3f', p_ks, D);
    lines{end+1} = sprintf('Probability of superiority PS = %.3f (bootstrap 95%% CI: [%.3f, %.3f])', PS, PS_CI(1), PS_CI(2));
    lines{end+1} = sprintf('Median(Q): in = %.3f | out = %.3f | Δ = %.3f', med_in, med_out, dmed);
    lines{end+1} = sprintf('Frac(Q >= %.3f): in = %.3f | out = %.3f', wet_thr0, frac_in, frac_out);
    if isfield(robustness,'did_run') && robustness.did_run
        lines{end+1} = sprintf('Robustness (stride=%d): ranksum p=%.3g | KS p=%.3g | PS_rank=%.3f', ...
            robustness.stride_k, robustness.p_mw, robustness.p_ks, robustness.PS_rank);
    else
        lines{end+1} = sprintf('Robustness: not run or insufficient subsample');
    end

    text(0.02, 0.95, lines, 'VerticalAlignment','top', 'FontName','Courier', 'FontSize', 11);
    exportgraphics(fig, fullfile(pngdir,'stats_summary_PIXELS.png'), 'Resolution', 220);
    close(fig);

end

function M = confusion_metrics_vector(wet, sink)
% wet, sink: logical vectors on same domain
    wet = logical(wet);
    sink = logical(sink);

    TP = nnz(wet & sink);
    FP = nnz(wet & ~sink);
    FN = nnz(~wet & sink);
    TN = nnz(~wet & ~sink);

    TPR = TP / max(TP+FN, 1);
    FPR = FP / max(FP+TN, 1);
    precision = TP / max(TP+FP, 1);
    specificity = TN / max(TN+FP, 1);
    balAcc = 0.5*(TPR + specificity);
    F1 = 2*TP / max(2*TP + FP + FN, 1);

    denom = sqrt( max(double((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)), 1) );
    MCC = (TP*TN - FP*FN) / denom;

    OR = ((TP+0.5)*(TN+0.5)) / ((FP+0.5)*(FN+0.5));

    M = struct('TP',TP,'FP',FP,'FN',FN,'TN',TN, ...
               'TPR',TPR,'FPR',FPR,'precision',precision,'specificity',specificity, ...
               'balAcc',balAcc,'MCC',MCC,'F1',F1,'oddsRatio',OR);
end

function [Qsub, sinksub] = subsample_stride(Q, sink_mask, valid, stride_k)
% Spatial subsampling to reduce autocorrelation:
% Keep every stride_k pixel in both row/col directions, restricted to valid.
%
% Returns vectors Qsub and logical sinksub on that subsampled set.

    if stride_k <= 1
        idx = find(valid);
        Qsub = Q(idx);
        sinksub = sink_mask(idx);
        return
    end

    [nr,nc] = size(Q);
    rr = 1:stride_k:nr;
    cc = 1:stride_k:nc;

    V = false(nr,nc);
    V(rr,cc) = true;
    V = V & valid;

    idx = find(V);
    Qsub = Q(idx);
    sinksub = sink_mask(idx);
end

function [PS, CI, psb] = prob_superiority_mc(qIn, qOut, cfg)
% PS = P(qIn > qOut) + 0.5*P(qIn == qOut)
% Monte Carlo estimate + bootstrap CI using paired random sampling.
%
% Usage:
%   [PS, CI] = prob_superiority_mc(qIn, qOut, cfg)
%
% cfg fields (all optional):
%   nPairs      : number of pairs for point estimate (default 2e6)
%   B           : bootstrap replicates (default 300)
%   nPairsBoot  : pairs per bootstrap replicate (default 3e5)
%   seed        : RNG seed (default 1)
%
% Outputs:
%   PS  : point estimate
%   CI  : [2.5 97.5] percentile bootstrap CI
%   psb : bootstrap samples (B x 1)

    if nargin < 3 || isempty(cfg), cfg = struct(); end

    % defaults (version-agnostic)
    if ~isfield(cfg,'nPairs')     || isempty(cfg.nPairs),     cfg.nPairs = 2e6;   end
    if ~isfield(cfg,'B')         || isempty(cfg.B),          cfg.B = 300;         end
    if ~isfield(cfg,'nPairsBoot')|| isempty(cfg.nPairsBoot), cfg.nPairsBoot = 3e5;end
    if ~isfield(cfg,'seed')      || isempty(cfg.seed),       cfg.seed = 1;        end

    % basic checks
    qIn  = qIn(:);
    qOut = qOut(:);

    qIn  = qIn(isfinite(qIn));
    qOut = qOut(isfinite(qOut));

    nIn  = numel(qIn);
    nOut = numel(qOut);

    if nIn == 0 || nOut == 0
        error('prob_superiority_mc:EmptyInput', 'qIn and qOut must be nonempty after finite filtering.');
    end

    rng(cfg.seed);

    % Point estimate
    nPairs = min(double(cfg.nPairs), 1e9); % guard
    iIn  = randi(nIn,  nPairs, 1);
    iOut = randi(nOut, nPairs, 1);

    a = qIn(iIn);
    b = qOut(iOut);

    PS = mean(a > b) + 0.5*mean(a == b);

    % Bootstrap
    B = double(cfg.B);
    psb = nan(B,1);

    nPairsBoot = min(double(cfg.nPairsBoot), 1e9);
    for k = 1:B
        iIn  = randi(nIn,  nPairsBoot, 1);
        iOut = randi(nOut, nPairsBoot, 1);
        a = qIn(iIn);
        b = qOut(iOut);
        psb(k) = mean(a > b) + 0.5*mean(a == b);
    end

    CI = prctile(psb, [2.5 97.5]);
end
