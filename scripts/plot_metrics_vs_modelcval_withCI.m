function plot_metrics_vs_modelcval_withCI(S, cfg)
% G-mean vs model median WITH bootstrapped confidence intervals
% - Y = metric EXP from S.Gstats (default 'G_ISPBwet_OSPBdry_raw') with CI (and inner CIi)
% - X = model medians via resolve_model_medians_robust(S)
% - Flats faded; non-flats labeled; optional manual point shown
%
% Usage:
%   L = load('figs_out/artifacts/S_plot_small_YYYYMMDD_HHmmss.mat');
%   plot_metrics_vs_modelcval_withCI(L.S_plot, L.cfg_plot)

% ---------------- config ----------------
if nargin < 2 || isempty(cfg), cfg = struct(); end
cfg = setdef(cfg, 'outdir', '.');
cfg = setdef(cfg, 'run_id', char(datetime('now','Format','yyyyMMdd_HHmmss')));
if ~isfield(cfg,'region') || ~isfield(cfg.region,'mode') || isempty(cfg.region.mode)
    cfg.region.mode = 'ALL';
end

% Default metric is the new regional G-mean
% 'G_ISPBwet_OSPBdry_raw' | 'REC_ISPB_wet_raw' | 'REC_OSPB_dry_raw'
% plus legacy options: 'G_raw' | 'REC_raw' | 'SPEC_raw' | 'ACC_raw' | 'PREC_raw' | 'F1_raw'
cfg = setdef(cfg, 'metric', 'G_ISPBwet_OSPBdry_raw');

cfg = setdef(cfg, 'label_nonflats', true);
cfg = setdef(cfg, 'exclude_names', {'MeanNonFlat'});
cfg = setdef(cfg, 'draw_inner_band', false);   % draw 50% CIi as thicker whisker
cfg = setdef(cfg, 'manual_point', struct( ...
    'enable', false, ...
    'x', 46.59, ...
    'y', 0.846, ...
    'label', 'MCMC', ...
    ... % Optional CIs you can manually type in (should match chosen metric)
    'y_CI95',   [NaN NaN], ...   % e.g., [0.84 0.91]
    'y_CI50',   [NaN NaN], ...       % optional inner/"thick" band
    'x_CI95',   [NaN NaN], ...       % optional horizontal CI on x
    'x_CI50',   [NaN NaN], ...       % optional inner x band
    'capsize',  2, ...
    'linew',    1.2, ...
    'band',     true, ...
    'alpha',    0.22 ));

cfg = setdef(cfg, 'fonts', struct('axis',11,'label',12));
cfg = setdef(cfg, 'colors', struct('flat',[0.7 0.7 0.7], 'non',[0.00 0.45 0.70], 'manual',[0.90 0.20 0.10]));

% ---------------- pull metric (+ CIs) ----------------
[Y, Ylo, Yhi, Ylo_i, Yhi_i, ok] = pull_metric_with_ci(S, cfg.metric);
assert(ok, 'Metric "%s" with CI not found in S.Gstats. Did bootstrap run?', cfg.metric);

% ---------------- X axis (model medians; robust fallbacks) ------------
X = resolve_model_medians_robust(S);

% ---------------- names / flats / excludes ----------------------------
nM = numel(X);
names   = repmat("", nM, 1);
titles  = names;
is_flat = false(nM,1);

if isfield(S,'names') && ~isempty(S.names)
    names = string(S.names(:));
    is_flat = startsWith(names, 'Flat_');
end
if isfield(S,'titles') && ~isempty(S.titles) && numel(S.titles)==nM
    titles = string(S.titles(:));
else
    titles = names;
end

% exclude specific aggregate names (e.g., 'MeanNonFlat')
excl = false(size(names));
if ~isempty(cfg.exclude_names)
    lowNames = lower(strtrim(names));
    lowEx    = lower(strtrim(string(cfg.exclude_names(:))));
    for i=1:numel(lowEx)
        excl = excl | strcmp(lowNames, lowEx(i));
    end
end

% order, masks
[x_ord, ord] = sort(X(:));
Y_ord     = Y(ord);
Ylo_ord   = Ylo(ord);   Yhi_ord   = Yhi(ord);
Yloi_ord  = Ylo_i(ord); Yhii_ord  = Yhi_i(ord);
names_ord = names(ord);
titles_ord= titles(ord);
isflat_ord= is_flat(ord);
excl_ord  = excl(ord);

good = isfinite(x_ord) & isfinite(Y_ord) & ~excl_ord;
xspan = max(1, range_nonempty(x_ord(good)));

% ---------------- plot ----------------
f = figure('Color','w','Name','Metric vs Model Median'); 
ax = axes('Parent',f); hold(ax,'on'); box(ax,'off');

% Flat models (faded)
maskF = isflat_ord & good;

% Non-flat models (colored)
maskN = ~isflat_ord & good;
if any(maskN)
    draw_ci_errorbars(ax, x_ord(maskN), Y_ord(maskN), Ylo_ord(maskN), Yhi_ord(maskN), ...
        'Color', cfg.colors.non, 'CapSize', 2, 'LineW', 1.0, 'Band', true, 'Alpha', 0.25, ...
        'XSpan', xspan);

    if cfg.draw_inner_band && all(isfinite(Yloi_ord(maskN)))
        draw_ci_errorbars(ax, x_ord(maskN), Y_ord(maskN), Yloi_ord(maskN), Yhii_ord(maskN), ...
            'Color', cfg.colors.non, 'CapSize', 0, 'LineW', 1.6, 'Band', false);
    end
    hPub = scatter(ax, x_ord(maskN), Y_ord(maskN), 42, 'o', ...
        'MarkerFaceColor',cfg.colors.non,'MarkerEdgeColor','k', ...
        'LineWidth',0.7,'DisplayName','Published/geophysical GHF models');
end

% ---- Interpolated trend for flat models (dashed line) ----
if any(maskF)
    [xq_flat, yq_flat] = build_flat_trend(x_ord(maskF), Y_ord(maskF));
    if ~isempty(xq_flat)
        hFlatTrend = plot(ax, xq_flat, yq_flat, '--', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 1.2, 'DisplayName', 'Uniform GHF trend');
    else
        hFlatTrend = [];
    end
else
    hFlatTrend = [];
end

% Optional manual point (with optional manual CIs)
hMan = [];
if isfield(cfg,'manual_point') && isstruct(cfg.manual_point) && isfield(cfg.manual_point,'enable') && cfg.manual_point.enable
    mp = cfg.manual_point;

    % --- Vertical (Y) CI bands for the manual point ---
    if isfield(mp,'y_CI95') && numel(mp.y_CI95)==2 && any(isfinite(mp.y_CI95))
        draw_ci_errorbars(ax, mp.x, mp.y, mp.y_CI95(1), mp.y_CI95(2), ...
            'Color', cfg.colors.manual, 'CapSize', mp.capsize, 'LineW', mp.linew, ...
            'Band', mp.band, 'Alpha', mp.alpha, 'XSpan', xspan);
    end
   
    if isfield(mp,'y_CI50') && numel(mp.y_CI50)==2 && all(isfinite(mp.y_CI50))
        draw_ci_errorbars(ax, mp.x, mp.y, mp.y_CI50(1), mp.y_CI50(2), ...
            'Color', cfg.colors.manual, 'CapSize', 0, 'LineW', mp.linew*1.4, ...
            'Band', false);
    end

    % --- Manual point marker + label ---
    hMan = scatter(ax, mp.x, mp.y, 52, 'o', 'MarkerFaceColor',cfg.colors.manual, ...
        'MarkerEdgeColor','k', 'LineWidth',1.0, 'DisplayName', mp.label);

    dx = 0.01 * max(1, range_nonempty(x_ord(good)));
    text(ax, mp.x + dx, mp.y, ['  ' mp.label], ...
        'Interpreter','none','FontSize',10, 'FontWeight','bold', ...
        'Color',[0.60 0 0], 'VerticalAlignment','bottom');
end

% Labels for non-flats
if cfg.label_nonflats && any(maskN)
    dx = 0.01 * max(1, range_nonempty(x_ord(good)));
    nn = titles_ord(maskN); xn = x_ord(maskN); yn = Y_ord(maskN);
    for k = 1:numel(xn)
        if strlength(nn(k)) > 0
            text(ax, xn(k)+dx, yn(k), char(nn(k)), ...
                'Interpreter','none','FontSize',12, ...
                'VerticalAlignment','middle','Color',[0.05 0.05 0.05]);
        end
    end
end

xlabel(ax, 'Median geothermal heat flow (mW m^{-2})','Interpreter','tex','FontSize',cfg.fonts.label);
ylabel(ax, metric_label(cfg.metric), 'Interpreter','tex','FontSize',cfg.fonts.label);
ylim(ax,[0 1]); grid(ax,'on'); grid(ax,'minor'); set(ax,'FontSize',cfg.fonts.axis);

% Robust legend building (only include handles that exist)
Lh = []; Ltxt = {};
if ~isempty(hMan) && isgraphics(hMan)
    Lh(end+1)   = hMan; 
    Ltxt{end+1} = cfg.manual_point.label; 
end
if exist('hPub','var') && ~isempty(hPub) && isgraphics(hPub)
    Lh(end+1)   = hPub;
    Ltxt{end+1} = 'Published/geophysical models';
end
if exist('hFlat','var') && ~isempty(hFlat) && isgraphics(hFlat)
    Lh(end+1)   = hFlat;
    Ltxt{end+1} = 'Uniform GHF fields';
end
if exist('hFlatTrend','var') && ~isempty(hFlatTrend) && isgraphics(hFlatTrend)
    Lh(end+1)   = hFlatTrend;
    Ltxt{end+1} = 'Uniform GHF trend';
end
if ~isempty(Lh)
    legend(ax, Lh, Ltxt, 'Location','northeast','Box','on','Orientation','vertical');
end

xlim([20 80]);
base = fullfile(cfg.outdir, sprintf('%s_vs_ModelMedian_CI_%s_%s', cfg.metric, upper(cfg.region.mode), cfg.run_id));
exportgraphics(f, [base '.png'], 'Resolution',300);
fprintf('[plot] saved: %s.png\n', base);

end

% ================= helpers =================
function cfg = setdef(cfg, field, val)
    if ~isfield(cfg,field) || isempty(cfg.(field)), cfg.(field) = val; end
end

function [Y, lo, hi, loi, hii, ok] = pull_metric_with_ci(S, field)
% field like 'G_raw', 'REC_raw', 'G_ISPBwet_OSPBdry_raw', etc.
% Expects .EXP, .CI, .CIi under S.Gstats.(field)
Y = []; lo = []; hi = []; loi = []; hii = []; ok = false;
if ~isfield(S,'Gstats') || ~isfield(S.Gstats, field), return; end
Gf = S.Gstats.(field);
if ~isfield(Gf,'EXP') || ~isfield(Gf,'CI') || isempty(Gf.CI), return; end
Y  = Gf.EXP(:);
CI = Gf.CI;
lo = CI(:,1); hi = CI(:,2);
if isfield(Gf,'CIi') && ~isempty(Gf.CIi)
    CIi = Gf.CIi; loi = CIi(:,1); hii = CIi(:,2);
else
    loi = NaN(size(lo)); hii = NaN(size(hi));
end
ok = true;
end

function X = resolve_model_medians_robust(S)
% Preferred order:
% 1) S.model_cvals_median
% 2) S.model_quantiles.q50
% 3) S.model_cvals (means)
% 4) compute from S.models or S.models_ds (if present)

if isfield(S,'model_cvals_median') && ~isempty(S.model_cvals_median)
    X = S.model_cvals_median(:); return;
end
if isfield(S,'model_quantiles') && isstruct(S.model_quantiles) ...
        && isfield(S.model_quantiles,'q50') && ~isempty(S.model_quantiles.q50)
    X = S.model_quantiles.q50(:); return;
end
if isfield(S,'model_cvals') && ~isempty(S.model_cvals)
    X = S.model_cvals(:); return;
end

% compute from models (full) or downsampled cache
flds = {};
src  = '';
if isfield(S,'models') && isstruct(S.models)
    flds = fieldnames(S.models); src = 'models';
elseif isfield(S,'models_ds') && isstruct(S.models_ds)
    flds = fieldnames(S.models_ds); src = 'models_ds';
else
    error('resolve_model_medians_robust:NeedModelCentralValues: Provide S.model_cvals(_median) or include S.models/_ds.');
end
nM = numel(flds); X = nan(nM,1);
roi_idx = [];
if isfield(S,'roi_mask') && ~isempty(S.roi_mask), roi_idx = find(S.roi_mask); end
for i = 1:nM
    M = S.(src).(flds{i});
    vals = double(M(:));
    if ~isempty(roi_idx) && isequal(numel(vals), numel(S.roi_mask))
        vals = double(M(roi_idx));
    end
    X(i) = median(vals, 'omitnan');
end
end

function draw_ci_errorbars(ax, x, y, ylo, yhi, varargin)
% Vertical CIs (whiskers on Y) at X positions.

p = inputParser;
addParameter(p,'Color',[0 0 0],@isnumeric);
addParameter(p,'CapSize',6,@isnumeric);
addParameter(p,'LineW',1.0,@isnumeric);
addParameter(p,'Band',false,@islogical);
addParameter(p,'Alpha',0.25,@isnumeric);
addParameter(p,'XSpan',[],@(v)isnumeric(v)&&isscalar(v)&&isfinite(v)&&v>0);  % NEW
parse(p, varargin{:});
C  = p.Results.Color; cs = p.Results.CapSize; lw = p.Results.LineW;
useBand = p.Results.Band; A = p.Results.Alpha; XSpan = p.Results.XSpan;

% Use provided XSpan if given; else try XLim; else fall back to data span
if ~isempty(XSpan)
    xr = XSpan;
else
    xl = get(ax,'XLim'); xr = diff(xl);
    if ~(isfinite(xr) && xr>0)
        xr = max(x(:)) - min(x(:));
        if ~(isfinite(xr) && xr>0), xr = 1; end
    end
end

capHalf  = (cs/200) * xr;
bandHalf = 0.0025 * xr;

for k = 1:numel(x)
    if ~(isfinite(x(k)) && isfinite(y(k)) && isfinite(ylo(k)) && isfinite(yhi(k))), continue; end
    if ylo(k) > yhi(k), [ylo(k), yhi(k)] = deal(yhi(k), ylo(k)); end

    % band first so whiskers/caps sit on top
    if useBand
        xv = [x(k)-bandHalf, x(k)+bandHalf, x(k)+bandHalf, x(k)-bandHalf];
        yv = [ylo(k),        ylo(k),        yhi(k),        yhi(k)];
        patch('Parent',ax,'XData',xv,'YData',yv, ...
              'FaceColor',C,'FaceAlpha',A,'EdgeColor','none','HandleVisibility','off');
    end

    line(ax, [x(k) x(k)], [ylo(k) yhi(k)], 'Color', C, 'LineWidth', lw, ...
         'Clipping','off','HandleVisibility','off');

    if cs > 0
        line(ax, [x(k)-capHalf x(k)+capHalf], [ylo(k) ylo(k)], 'Color', C, 'LineWidth', lw, ...
             'Clipping','off','HandleVisibility','off');
        line(ax, [x(k)-capHalf x(k)+capHalf], [yhi(k) yhi(k)], 'Color', C, 'LineWidth', lw, ...
             'Clipping','off','HandleVisibility','off');
    end
end

set(ax,'Layer','top');
end

function r = range_nonempty(v)
v = v(:); v = v(isfinite(v));
r = 1; if ~isempty(v), r = max(v)-min(v); end
end

function s = metric_label(m)
m = string(m);
switch m
    case "G_raw"
        s = 'G-mean (global REC×SPEC, bootstrapped)';
    case "REC_raw"
        s = 'Recall (all sinks, bootstrapped)';
    case "SPEC_raw"
        s = 'Specificity (all sinks, bootstrapped)';
    case "ACC_raw"
        s = 'Accuracy (all sinks, bootstrapped)';
    case "PREC_raw"
        s = 'Precision (all sinks, bootstrapped)';
    case "F1_raw"
        s = 'F1 score (all sinks, bootstrapped)';

    % ---- NEW regional metrics ----
    case "REC_ISPB_wet_raw"
        s = 'Recall of wet sinks in ISPB (bootstrapped)';
    case "REC_OSPB_dry_raw"
        s = 'Recall of dry sinks in OSPB (bootstrapped)';
    case "G_ISPBwet_OSPBdry_raw"
        s = 'G-mean: ISPB wet recall × OSPB dry recall (bootstrapped)';

    otherwise
        s = char(m);
end
end

function [xq, yq] = build_flat_trend(x, y)
% Build a smooth, monotone-in-x trend for flat models
% - Sorts & uniques x
% - Applies a gentle moving-median denoise if enough points
% - Interpolates with 'pchip' to a dense grid for a smooth dashed line
    xq = []; yq = [];
    ok = isfinite(x) & isfinite(y);
    x  = x(ok); y = y(ok);
    if numel(x) < 2, return; end

    [xs, idx] = sort(x); ys = y(idx);
    [xu, ia]  = unique(xs, 'stable'); yu = ys(ia);

    % gentle smoothing only if lots of flats; otherwise just use yu
    if numel(xu) >= 5
        w = max(3, round(numel(xu)/8));                % gentle smoothing
        yu_s = movmedian(yu, w, 'omitnan');
    else
        yu_s = yu;
    end

    xq = linspace(min(xu), max(xu), 300);
    yq = interp1(xu, yu_s, xq, 'pchip', 'extrap');
end

function draw_hci_errorbars(ax, x, y, xlo, xhi, varargin)
% Horizontal CIs (whiskers on X) at Y position.
p = inputParser;
addParameter(p,'Color',[0 0 0],@isnumeric);
addParameter(p,'CapSize',6,@isnumeric);
addParameter(p,'LineW',1.0,@isnumeric);
addParameter(p,'Band',false,@islogical);
addParameter(p,'Alpha',0.25,@isnumeric);
parse(p, varargin{:});
C = p.Results.Color; cs = p.Results.CapSize; lw = p.Results.LineW;
useBand = p.Results.Band; A = p.Results.Alpha;

% derive small widths in data units (for cap height)
yr = diff(get(ax,'YLim')); if yr<=0 || ~isfinite(yr), yr = 1; end
capHalf = (cs/200) * yr;       % cap half-height in data units
bandHalf= 0.0025 * yr;         % thin band thickness

if ~(isfinite(x) && isfinite(y) && isfinite(xlo) && isfinite(xhi)), return; end
if xlo > xhi, tmp=xlo; xlo=xhi; xhi=tmp; end

% main horizontal whisker
line(ax, [xlo xhi], [y y], 'Color', C, 'LineWidth', lw);

% caps
if cs > 0
    line(ax, [xlo xlo], [y-capHalf y+capHalf], 'Color', C, 'LineWidth', lw);
    line(ax, [xhi xhi], [y-capHalf y+capHalf], 'Color', C, 'LineWidth', lw);
end

% optional skinny band
if useBand
    xv = [xlo, xhi, xhi, xlo];
    yv = [y-bandHalf, y-bandHalf, y+bandHalf, y+bandHalf];
    patch('Parent',ax,'XData',xv,'YData',yv,'FaceColor',C,'FaceAlpha',A,'EdgeColor','none');
end
end
