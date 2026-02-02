%% ================================================================
%  GHF model evaluation (v3) — weighted specularity truth + full outputs
%  - Preserves legacy binary point-metric outputs (hard truth)
%  - Adds weighted-truth point metrics (soft truth)
%  - Runs analytic expected metrics + bootstrap on cfg.truth.mode ("binary" or "weighted")
%  - Writes CI summary CSV (like original)
%  - Saves S_plot_small bundle (like original, plus y_bin/y_w fields)
%  - Caches per-model full-grid GΔ fields (like original)
%
%  Drop-in replacement for: bootstrap_ghf_component.m
%  Save as: bootstrap_ghf_component_weighted.m
%  ================================================================
function S = bootstrap_ghf_component_weighted(cfg_override)

check_dependencies(struct('autoInstall',true,'needAMT',false,'verbose',true));
addpath('datasets_for_gmin');
dbstop if error
run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));

%% -------------------- CONFIG --------------------
cfg = struct( ...
  'spec_thresh',          0.2, ...        % legacy binary threshold (kept)
  'v_keep',               10,  ...        % ice velocity threshold (m/yr)
  'decision_margin',      0,   ...        % wet if (GΔ >= margin)
  'decision_margin_mode', 'fixed', ...    % fixed | auto_mean | auto_median | auto_per_model_mean | auto_per_model_median
  'outdir',               'figs_out', ...
  'overwrite',            true, ...
  'overlay_alpha',        0.4, ...
  'font_size',            20, ...
  'to_single',            true, ...
  'sweep_mode',           true ...
);

% Truth label mode (binary or weighted)
cfg.truth = struct( ...
  'mode', 'weighted', ...              % 'binary' | 'weighted'
  'weight_mode', 'logistic', ...       % 'ramp' | 'logistic' | 'power'
  ... % calibration targets for logistic auto-calibration
  'auto', struct( ...
      'enable', true, ...
      'q_low',  0.60, ...   % e.g., median of ROI Sc
      'q_high', 0.95, ...   % strong tail of ROI Sc
      'w_low',  0.60, ...   % desired weight at Sc_low
      'w_high', 0.95  ...   % desired weight at Sc_high
  ), ...
  ... % optional manual params (used if auto disabled)
  't0',   0.2, ...
  'tmid', 0.6, ...
  'k',    12, ...
  'p',    2 ...
);

cfg.save_artifacts = true;
cfg.artifact_dir   = fullfile(cfg.outdir,'artifacts');
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end
cfg.run_id = run_id;

% Region split
cfg.region = struct( ...
  'mode',          'ALL', ...         % 'ISPB' | 'OSPB' | 'SPB' | 'ALL'
  'left_label',    'ISPB', ...
  'right_label',   'OSPB', ...
  'cache',         true ...
);

% Sinks
cfg.sinks = struct('marker','.', 'size',4, 'color',[0.05 0.05 0.05], 'alpha',0.85);

cfg.skip_static_build = true;

% Uncertainty config (same spirit as original)
cfg.uncertainty = struct( ...
  'mode','analytic', ...
  'bootstrap_mode','pixel', ...    % 'pixel' | 'component'
  'n_boot',1000, ...
  'seed',42, ...
  'mfrac',0.7, ...
  'band_main',0.95, ...
  'band_inner',0.50, ...
  'noise_enable',false, ...
  'noise_K',1, ...
  'pred_mode','bernoulli', ...     % 'bernoulli' | 'noisy_threshold'
  'class_balance',false, ...       % only meaningful for HARD labels
  'pps_components',false ...       % component bootstrap: PPS by component size
);

cfg.plot_labels = struct('avoid_overlap',true,'max_iter',200,'step',0.01,'pad',0.005);
cfg.synthetic_only = false;

if nargin > 0 && ~isempty(cfg_override)
    cfg = merge_structs(cfg, cfg_override);
end
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end

%% =========================== Load data ===========================
fprintf('Loading data ...\n');
datafile = fullfile('datasets_for_gmin','gmin_data.mat');
if exist(datafile,'file')
    fprintf('[bootstrap] Loading %s ... ', datafile);
    L = load(datafile,'S'); S = L.S; clear L; fprintf('ok.\n');
else
    warning('Data file not found. Building now...');
    S = build_data_core('to_single',true,'save',true, ...
                    'outdir','datasets_for_gmin','outfile','gmin_data.mat');
    S = add_g_fields(fullfile('datasets_for_gmin','gmin_data.mat'), ...
                 'advection_file','datasets_for_gmin/longitudinal_advection_maps.mat', ...
                 'advection_var','Delta_q_adv', 'advection_units','mWm2', ...
                 'save',true, 'outfile', fullfile('datasets_for_gmin','gmin_data.mat'));
    S = add_g_uncertainty('datasets_for_gmin/gmin_data.mat', ...
                'sigma_Ts',1.5,'sigma_Mb',0.02,'sigma_H',50, ...
                'delta_Ts',0.5,'delta_Mb',0.01,'delta_H',5, ...
                'cache_dir','', 'use_cache',true, 'save',true,'outfile','');
end

% Light safety fields
if ~isfield(S,'X_km') || ~isfield(S,'Y_km')
    S.X_km = S.Xgrid/1000; S.Y_km = S.Ygrid/1000;
end
if ~isfield(S,'spec_invalid') && isfield(S,'Q')
    S.spec_invalid = (S.Q == 9999);
end
fprintf('[bootstrap] Grid: %dx%d | models: %d\n', size(S.Xgrid,1), size(S.Xgrid,2), numel(S.names));
nM = numel(S.names);
S.model_cvals = nan(nM,1);

%% =========================== Load SPB masks ===========================
[nr,nc] = size(S.Xgrid);
REG_MASK = true(nr,nc);

mode = upper(string(cfg.region.mode));
maskfile = fullfile('datasets_for_gmin','spb_masks.mat');

if exist(maskfile,'file') ~= 2
    warning('[REG_MASK] %s missing -> using ALL (full grid).', maskfile);
    S.rect_mask = true(nr,nc);
else
    m = load(maskfile);

    mask_ISPB = [];
    if isfield(m,'mask_ISPB_c'), mask_ISPB = full(m.mask_ISPB_c);
    elseif isfield(m,'mask_ISPB'), mask_ISPB = m.mask_ISPB; end

    mask_OSPB = [];
    if isfield(m,'mask_OSPB_c'), mask_OSPB = full(m.mask_OSPB_c);
    elseif isfield(m,'mask_OSPB'), mask_OSPB = m.mask_OSPB; end

    mask_SPB  = [];
    if isfield(m,'mask_SPB_c'), mask_SPB = full(m.mask_SPB_c);
    elseif isfield(m,'mask_SPB'), mask_SPB = m.mask_SPB; end

    rect = true(nr,nc);
    if isfield(m,'rect_mask'), rect = m.rect_mask; end

    S.mask_ISPB = expand_or_resize_mask(mask_ISPB, m, size(S.Xgrid), 'mask_ISPB');
    S.mask_OSPB = expand_or_resize_mask(mask_OSPB, m, size(S.Xgrid), 'mask_OSPB');
    S.mask_SPB  = expand_or_resize_mask(mask_SPB,  m, size(S.Xgrid), 'mask_SPB');
    S.rect_mask = expand_or_resize_mask(rect,      m, size(S.Xgrid), 'rect_mask');

    switch mode
        case "ISPB", REG_MASK = S.mask_ISPB & S.rect_mask;
        case "OSPB", REG_MASK = S.mask_OSPB & S.rect_mask;
        case "SPB",  REG_MASK = S.mask_SPB  & S.rect_mask;
        case "ALL",  REG_MASK = S.rect_mask;
        otherwise, error('cfg.region.mode must be ISPB | OSPB | SPB | ALL');
    end
end

if ~isequal(size(REG_MASK), size(S.Xgrid))
    warning('[REG_MASK] size mismatch; resizing nearest.');
    REG_MASK = logical(imresize(REG_MASK, size(S.Xgrid), 'nearest'));
end
fprintf('[REG_MASK] mode=%s | nnz=%d\n', mode, nnz(REG_MASK));

%% ========== ROI masks (sinks + slow flow + region) ==========
if exist('datasets_for_gmin/sink_mask_new.mat','file')==2 && exist('datasets_for_gmin/sink_mask_comp.mat','file')==2
    t = load('datasets_for_gmin/sink_mask_new.mat', 'sink_mask'); S.sink_mask = t.sink_mask; clear t;
    t = load('datasets_for_gmin/sink_mask_comp.mat', 'comp_id','sink_mask_comp');
    S.comp_id = t.comp_id; S.sink_mask_comp = t.sink_mask_comp; clear t;
else
    make_sink_mask();
    t = load('datasets_for_gmin/sink_mask_new.mat', 'sink_mask'); S.sink_mask = t.sink_mask; clear t;
    t = load('datasets_for_gmin/sink_mask_comp.mat', 'comp_id','sink_mask_comp');
    S.comp_id = t.comp_id; S.sink_mask_comp = t.sink_mask_comp; clear t;
end

valid_mask = isfinite(S.Q) & ~S.spec_invalid & isfinite(S.H) & isfinite(S.icevel);
slow_mask  = valid_mask & ~(S.icevel > cfg.v_keep);

S.viz_mask = slow_mask & REG_MASK;
S.roi_mask = slow_mask & REG_MASK & S.sink_mask_comp;

idx     = find(S.roi_mask);
idx_med = find(S.viz_mask);  % median domain: slow & REG, no sinks

S.eval_mask_full = S.roi_mask;
S.eval_mask      = true(numel(idx),1);

% Regional vectors
mask_ISPB_vec = false(numel(idx),1);
mask_OSPB_vec = false(numel(idx),1);
if isfield(S,'mask_ISPB'), mask_ISPB_vec = S.mask_ISPB(idx); end
if isfield(S,'mask_OSPB'), mask_OSPB_vec = S.mask_OSPB(idx); end
S.mask_ISPB_vec = mask_ISPB_vec;
S.mask_OSPB_vec = mask_OSPB_vec;

%% ==================== S.Q statistics (computed here) ====================
S.spec_stats = struct();
S.spec_stats.region_mode = char(mode);
S.spec_stats.v_keep = cfg.v_keep;

mask_roi = S.roi_mask;
mask_viz = S.viz_mask;
mask_reg = REG_MASK & valid_mask;

S.spec_stats.ROI = summarize_SQ(S.Q, mask_roi);
S.spec_stats.VIZ = summarize_SQ(S.Q, mask_viz);
S.spec_stats.REG = summarize_SQ(S.Q, mask_reg);

fprintf('\n[S.Q stats] ROI: N=%d | p50=%.3f p90=%.3f p95=%.3f p99=%.3f | mean=%.3f std=%.3f\n', ...
    S.spec_stats.ROI.N, S.spec_stats.ROI.p50, S.spec_stats.ROI.p90, S.spec_stats.ROI.p95, S.spec_stats.ROI.p99, ...
    S.spec_stats.ROI.mean, S.spec_stats.ROI.std);

%% ==================== Truth labels: binary + weighted ====================
Sc_roi = S.Q(idx);

% legacy binary truth (hard)
y_bin = double(isfinite(Sc_roi) & (Sc_roi ~= 9999) & (Sc_roi > cfg.spec_thresh));

% ---- Auto-calibrate logistic parameters from ROI anchors ----
cfg_used = cfg;

Sc = double(Sc_roi);
Sc = Sc(isfinite(Sc) & Sc~=9999);

if isfield(cfg.truth,'auto') && isfield(cfg.truth.auto,'enable') && cfg.truth.auto.enable
    a = cfg.truth.auto;

    Sc_low  = quantile(Sc, a.q_low);
    Sc_high = quantile(Sc, a.q_high);

    w_low  = max(min(a.w_low,  0.999), 0.001);
    w_high = max(min(a.w_high, 0.999), 0.001);

    L1 = log(w_low  /(1-w_low));
    L2 = log(w_high /(1-w_high));

    k_auto    = (L2 - L1) / max(Sc_high - Sc_low, eps);
    tmid_auto = Sc_low - L1 / max(k_auto, eps);

    cfg_used.truth.weight_mode = 'logistic';
    cfg_used.truth.k    = k_auto;
    cfg_used.truth.tmid = tmid_auto;

    S.truth_cal = struct('Sc_low',Sc_low,'Sc_high',Sc_high,'w_low',w_low,'w_high',w_high, ...
                         'tmid',tmid_auto,'k',k_auto, 'q_low',a.q_low,'q_high',a.q_high);
end

% weighted truth (soft)
y_w = spec_weight(Sc_roi, cfg_used.truth);

S.y_bin_vec = y_bin(:);
S.y_w_vec   = y_w(:);

% store full grids too (useful for plots)
S.y_bin_full = nan(size(S.Q), 'single');
S.y_w_full   = nan(size(S.Q), 'single');
S.y_bin_full(idx) = single(S.y_bin_vec);
S.y_w_full(idx)   = single(S.y_w_vec);

% pick which truth to USE downstream for analytic+bootstrap + plotting bundle
truth_mode = lower(string(cfg.truth.mode));
switch truth_mode
    case "weighted"
        y_use = S.y_w_vec;
    case "binary"
        y_use = S.y_bin_vec;
    otherwise
        error('cfg.truth.mode must be "binary" or "weighted"');
end
S.truth_mode_used = char(truth_mode);
S.y_use_vec = y_use(:);

% Mirror legacy naming expected by original downstream blocks
S.y_raw_vec  = S.y_use_vec;            % ROI-vector labels (hard or soft)
S.y_raw_full = S.y_w_full;             % full-grid label for plotting
if strcmpi(S.truth_mode_used,'binary')
    S.y_raw_full = S.y_bin_full;
end

%% ==================== Component grouping (for component bootstrap) ====================
S.boot_groups = struct('nC',0,'groups',{{}},'uniqC',[]);
if strcmpi(cfg.uncertainty.bootstrap_mode,'component')
    compM = S.comp_id(idx);
    if any(compM)
        [uniqC,~,gidx] = unique(compM(compM>0));
        groups = accumarray(gidx, find(compM>0), [], @(v){v});
        S.boot_groups.uniqC  = uniqC;
        S.boot_groups.groups = groups;
        S.boot_groups.nC     = numel(groups);
    end
end

%% ----------------- Per-model center values & GΔ cache -------------------
validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};
    Mi     = S.models.(fld);
    if startsWith(name_i,'Flat_')
        v = sscanf(name_i,'Flat_%f');
        S.model_cvals(i) = v;
    else
        S.model_cvals(i) = median(Mi(idx_med), 'omitnan');
    end
end

gdif_cache = cell(nM,1);
for i = 1:nM
    fld = validNames{i};
    gdif_cache{i} = S.Gdiff.(fld)(idx);  % ROI-vector
end

%% ==================== SANITY ====================
t0   = getfield_default(cfg_used.truth, 't0',   NaN);
tmid = getfield_default(cfg_used.truth, 'tmid', NaN);
k    = getfield_default(cfg_used.truth, 'k',    NaN);
p    = getfield_default(cfg_used.truth, 'p',    NaN);
wm   = getfield_default(cfg_used.truth, 'weight_mode', 'unknown');

fprintf('\n---------- SANITY ----------\n');
fprintf('[truth] mode_used=%s | spec_thresh(binarize)=%.2f | weight_mode=%s | t0=%.3f tmid=%.3f k=%.1f p=%.2f\n', ...
    S.truth_mode_used, cfg.spec_thresh, string(wm), t0, tmid, k, p);
fprintf('[truth] ROI base rate (binary) = %.3f%%\n', 100*mean(S.y_bin_vec>0));
fprintf('[truth] ROI mean weight        = %.3f\n', mean(S.y_w_vec, 'omitnan'));
fprintf('---------------------------\n\n');



%% ==================== POINT METRICS: legacy hard binary =================
fprintf('Evaluating *binary* point metrics on ROI (legacy hard truth)...\n');

TP  = zeros(nM,1); FP=TP; TN=TP; FN=TP;
ACC = TP; PREC=TP; REC=TP; F1=TP; SPEC=TP;
REC_ISPB_wet = nan(nM,1); REC_OSPB_dry = nan(nM,1); G_ISPBwet_OSPBdry = nan(nM,1);

for i = 1:nM
    x  = gdif_cache{i};
    Mv = S.eval_mask;
    dm_i = get_model_margin(cfg, S, i);
    prd  = x(Mv) >= dm_i;

    yr_bin = logical(S.y_bin_vec(Mv));
    m_raw  = compute_metrics_hard(prd, yr_bin);

    TP(i)=m_raw.TP; FP(i)=m_raw.FP; TN(i)=m_raw.TN; FN(i)=m_raw.FN;
    ACC(i)=m_raw.ACC; PREC(i)=m_raw.PREC; REC(i)=m_raw.REC; F1(i)=m_raw.F1; SPEC(i)=m_raw.SPEC;

    % regional binary recalls (legacy)
    pred_full = false(numel(idx),1); pred_full(Mv) = prd;
    wet_full  = logical(S.y_bin_vec);
    dry_full  = ~wet_full;

    mwet = mask_ISPB_vec & wet_full;
    mdry = mask_OSPB_vec & dry_full;

    if nnz(mwet)>0
        REC_ISPB_wet(i) = sum(pred_full(mwet)) / max(nnz(mwet), eps);
    end
    if nnz(mdry)>0
        REC_OSPB_dry(i) = sum(~pred_full(mdry)) / max(nnz(mdry), eps);
    end
    if isfinite(REC_ISPB_wet(i)) && isfinite(REC_OSPB_dry(i))
        G_ISPBwet_OSPBdry(i) = sqrt(max(REC_ISPB_wet(i),0) * max(REC_OSPB_dry(i),0));
    end
end

S.results_table_binary = table(string(S.names(:)), TP,FP,TN,FN, ACC,PREC,REC,F1,SPEC, ...
    REC_ISPB_wet, REC_OSPB_dry, G_ISPBwet_OSPBdry, ...
    'VariableNames', {'Model','TP','FP','TN','FN','ACC','PREC','REC','F1','SPEC', ...
                      'REC_ISPB_wet','REC_OSPB_dry','G_ISPBwet_OSPBdry'});

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
out_csv = fullfile(cfg.outdir, sprintf('results_binary_%s_%s.csv', char(mode), stamp));
writetable(S.results_table_binary, out_csv);
fprintf('Saved legacy binary results: %s\n', out_csv);

%% ==================== POINT METRICS: weighted truth (soft confusion) =====
fprintf('Evaluating *weighted-truth* point metrics on ROI (soft truth)...\n');

TPw  = zeros(nM,1); FPw=TPw; TNw=TPw; FNw=TPw;
ACCw = TPw; PRECw=TPw; RECw=TPw; F1w=TPw; SPECw=TPw;
RECw_ISPB = nan(nM,1); RECw_OSPB = nan(nM,1); Gw_reg = nan(nM,1);

yW_full = S.y_w_vec;
dryW_full = 1 - yW_full;

for i = 1:nM
    x  = gdif_cache{i};
    Mv = S.eval_mask;
    dm_i = get_model_margin(cfg, S, i);
    prd  = x(Mv) >= dm_i;

    yr_w = yW_full(Mv);
    m_w  = compute_metrics_soft(prd, yr_w);

    TPw(i)=m_w.TP; FPw(i)=m_w.FP; TNw(i)=m_w.TN; FNw(i)=m_w.FN;
    ACCw(i)=m_w.ACC; PRECw(i)=m_w.PREC; RECw(i)=m_w.REC; F1w(i)=m_w.F1; SPECw(i)=m_w.SPEC;

    % regional weighted recalls:
    pred_full = false(numel(idx),1); pred_full(Mv) = prd;

    mwet = mask_ISPB_vec;
    mdry = mask_OSPB_vec;

    denom_wet = sum(yW_full(mwet), 'omitnan');
    denom_dry = sum(dryW_full(mdry), 'omitnan');

    if denom_wet > 0
        RECw_ISPB(i) = sum(double(pred_full(mwet)).*yW_full(mwet), 'omitnan') / denom_wet;
    end
    if denom_dry > 0
        RECw_OSPB(i) = sum(double(~pred_full(mdry)).*dryW_full(mdry), 'omitnan') / denom_dry;
    end
    if isfinite(RECw_ISPB(i)) && isfinite(RECw_OSPB(i))
        Gw_reg(i) = sqrt(max(RECw_ISPB(i),0) * max(RECw_OSPB(i),0));
    end
end

S.results_table_weighted = table(string(S.names(:)), TPw,FPw,TNw,FNw, ACCw,PRECw,RECw,F1w,SPECw, ...
    RECw_ISPB, RECw_OSPB, Gw_reg, ...
    'VariableNames', {'Model','TPw','FPw','TNw','FNw','ACCw','PRECw','RECw','F1w','SPECw', ...
                      'RECw_ISPB_wet','RECw_OSPB_dry','Gw_ISPBwet_OSPBdry'});

out_csv = fullfile(cfg.outdir, sprintf('results_weighted_%s_%s.csv', char(mode), stamp));
writetable(S.results_table_weighted, out_csv);
fprintf('Saved weighted results: %s\n', out_csv);

%% ================= EXPECTED METRICS + BOOTSTRAP (truth-mode used) =========
fprintf('Analytic expected metrics + %s bootstrap on truth_mode_used="%s"...\n', ...
    cfg.uncertainty.bootstrap_mode, S.truth_mode_used);

stat_list = {'SPEC','REC','PREC','ACC','F1'};

% ROI-vector M-space
M        = S.eval_mask;                 % logical over ROI-vector
y_raw_M  = double(S.y_raw_vec(:));      % labels in ROI-vector
y_raw_M  = y_raw_M(M);
NM       = nnz(M);

% regional masks in M-space (for analytic + bootstrap)
mask_wet_ISPB_M = mask_ISPB_vec(M);   % subset mask; denominator handled below
mask_dry_OSPB_M = mask_OSPB_vec(M);

% σ(Gmin) per pixel in ROI-vector (0 if missing)
if isfield(S,'sigma_gmin') && ~isempty(S.sigma_gmin) && isequal(size(S.sigma_gmin), size(S.H))
    sigG_M_full = S.sigma_gmin(idx);
    sigG_M      = sigG_M_full(M);
else
    sigG_M = zeros(NM,1,'like',y_raw_M);
end

% Prepare σ(GΔ) aggregates inside ROI
S.sigma_gdif_mean   = nan(nM,1);
S.sigma_gdif_median = nan(nM,1);
scalar_sigma_per_model = nan(nM,1);

if isfield(cfg.uncertainty,'seed') && ~isempty(cfg.uncertainty.seed), rng(cfg.uncertainty.seed); end

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};

    sigM_grid = [];
    if isfield(S,'unc') && isfield(S.unc,fld) && isequal(size(S.unc.(fld)), size(S.H))
        sigM_grid = S.unc.(fld);
    elseif strcmp(name_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        sigM_grid = abs(S.bounds.Losing_max - S.bounds.Losing_min) ./ (2*1.96);
    end

    if ~isempty(sigM_grid)
        sigM_roi = sigM_grid(idx); sigM_roi = sigM_roi(M);
        sigma_gdel_vec = sqrt(sigM_roi.^2 + sigG_M.^2);
        S.sigma_gdif_mean(i)   = mean(sigma_gdel_vec, 'omitnan');
        S.sigma_gdif_median(i) = median(sigma_gdel_vec, 'omitnan');
        scalar_sigma_per_model(i) = mean(sigM_roi, 'omitnan');
    else
        guess = NaN;
        if isfield(S,'sigma_model_scalar') && numel(S.sigma_model_scalar)>=i
            guess = S.sigma_model_scalar(i);
        end
        if ~isfinite(guess), guess = 0; end
        scalar_sigma_per_model(i) = guess;

        S.sigma_gdif_mean(i)   = mean(sqrt((guess.^2) + (sigG_M.^2)), 'omitnan');
        S.sigma_gdif_median(i) = median(sqrt((guess.^2) + (sigG_M.^2)), 'omitnan');
    end
end

EXP_raw = struct(); CI_raw = struct(); CIi_raw = struct();
for s = stat_list
    EXP_raw.(s{1}) = nan(nM,1);
    CI_raw.(s{1})  = nan(nM,2);
    CIi_raw.(s{1}) = nan(nM,2);
end
EXP_raw.G = nan(nM,1); CI_raw.G = nan(nM,2); CIi_raw.G = nan(nM,2);

% Optional regional metrics for truth-mode used (works for hard or soft)
EXP_raw.REC_ISPB_wet      = nan(nM,1);
EXP_raw.REC_OSPB_dry      = nan(nM,1);
EXP_raw.G_ISPBwet_OSPBdry = nan(nM,1);
CI_raw.REC_ISPB_wet      = nan(nM,2);
CI_raw.REC_OSPB_dry      = nan(nM,2);
CI_raw.G_ISPBwet_OSPBdry = nan(nM,2);
CIi_raw.REC_ISPB_wet      = nan(nM,2);
CIi_raw.REC_OSPB_dry      = nan(nM,2);
CIi_raw.G_ISPBwet_OSPBdry = nan(nM,2);

use_comp_boot = strcmpi(cfg.uncertainty.bootstrap_mode,'component') && isfield(S,'boot_groups') && S.boot_groups.nC>0;
G_M = {}; if use_comp_boot, G_M = S.boot_groups.groups; end

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};
    gdif_M_all = gdif_cache{i}(M);

    % σ(model) vector for ROI-vector M-space
    sigM_M = zeros(NM,1,'like',gdif_M_all);
    if isfield(S,'unc') && isfield(S.unc, fld) && ~isempty(S.unc.(fld)) ...
            && isequal(size(S.unc.(fld)), size(S.H))
        tmp    = S.unc.(fld);
        sigM_M = tmp(idx); sigM_M = sigM_M(M);
    elseif strcmp(name_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        half_roi = 0.5*abs(S.bounds.Losing_max(idx) - S.bounds.Losing_min(idx));
        sigM_M   = (half_roi(M)) / 1.96;
    else
        sigM_M(:) = scalar_sigma_per_model(i);
    end

    dm_i  = get_model_margin(cfg, S, i);
    pw_all = wet_prob_with_margin(gdif_M_all, sigM_M, sigG_M, dm_i);
    V      = isfinite(pw_all) & isfinite(y_raw_M);
    pw = pw_all(V);
    yr = y_raw_M(V);
    Nf = numel(pw);
    if Nf==0, continue; end

    % regional subset masks within V
    mask_ISPB_V = mask_wet_ISPB_M(V);
    mask_OSPB_V = mask_dry_OSPB_M(V);

    % Analytic expectations (soft-label compatible)
    ETP  = sum(pw .* yr, 'omitnan');
    EFP  = sum(pw .* (1-yr), 'omitnan');
    EFN  = sum((1-pw) .* yr, 'omitnan');
    ETN  = sum((1-pw) .* (1-yr), 'omitnan');

    PREC_E = ETP / max(ETP + EFP, eps);
    SPEC_E = ETN / max(ETN + EFP, eps);
    REC_E  = ETP / max(ETP + EFN, eps);
    ACC_E  = (ETP + ETN) / max(ETP+EFP+EFN+ETN, eps);
    F1_E   = 2 * PREC_E * REC_E / max(PREC_E + REC_E, eps);

    EXP_raw.PREC(i)=PREC_E;
    EXP_raw.SPEC(i)=SPEC_E;
    EXP_raw.REC(i) =REC_E;
    EXP_raw.ACC(i) =ACC_E;
    EXP_raw.F1(i)  =F1_E;
    EXP_raw.G(i)   = sqrt(max(REC_E,0)*max(SPEC_E,0));

    % Regional expectations (soft-label compatible)
    % ISPB wet recall: numerator=sum(pw*y) denom=sum(y)
    denom_wet = sum(yr(mask_ISPB_V), 'omitnan');
    if denom_wet > 0
        EXP_raw.REC_ISPB_wet(i) = sum(pw(mask_ISPB_V).*yr(mask_ISPB_V), 'omitnan') / denom_wet;
    end
    % OSPB dry recall: numerator=sum((1-pw)*(1-y)) denom=sum(1-y)
    denom_dry = sum((1-yr(mask_OSPB_V)), 'omitnan');
    if denom_dry > 0
        EXP_raw.REC_OSPB_dry(i) = sum((1-pw(mask_OSPB_V)).*(1-yr(mask_OSPB_V)), 'omitnan') / denom_dry;
    end
    if isfinite(EXP_raw.REC_ISPB_wet(i)) && isfinite(EXP_raw.REC_OSPB_dry(i))
        EXP_raw.G_ISPBwet_OSPBdry(i) = sqrt(max(EXP_raw.REC_ISPB_wet(i),0) * max(EXP_raw.REC_OSPB_dry(i),0));
    end

    % ------------------- Bootstrap replicates -------------------
    B  = cfg.uncertainty.n_boot;
    BR = struct();
    for s = stat_list, BR.(s{1}) = zeros(B,1); end
    BR.G = zeros(B,1);
    BR.REC_ISPB_wet      = nan(B,1);
    BR.REC_OSPB_dry      = nan(B,1);
    BR.G_ISPBwet_OSPBdry = nan(B,1);

    pred_mode = lower(string(cfg.uncertainty.pred_mode));
    use_noisy = strcmp(pred_mode,"noisy_threshold") && isfield(cfg.uncertainty,'noise_enable') && cfg.uncertainty.noise_enable;
    s_all = []; % lazily compute

    if strcmpi(cfg.uncertainty.bootstrap_mode,'pixel')
        for b = 1:B
            % class_balance only for HARD labels
            uniqv = unique(yr(~isnan(yr)));
            isHard = all(ismember(uniqv, [0 1]));
            if cfg.uncertainty.class_balance && isHard
                pos = find(yr==1); neg = find(yr==0);
                nPos = numel(pos); nNeg = numel(neg);
                if nPos>0 && nNeg>0
                    rPos = pos(randi(nPos, nPos, 1));
                    rNeg = neg(randi(nNeg, nNeg, 1));
                    rb   = [rPos; rNeg];
                    rb   = rb(randperm(numel(rb)));
                else
                    rb = randi(Nf, Nf, 1);
                end
            else
                rb = randi(Nf, Nf, 1);
            end

            y_b = yr(rb);

            if use_noisy
                if isempty(s_all), s_all = sqrt(max(sigM_M.^2 + sigG_M.^2, eps)); end
                epsi = randn(numel(rb),1) .* (cfg.uncertainty.noise_K .* s_all(V));
                gjit = gdif_M_all; gjit = gjit(V);
                gjit = gjit(rb) + epsi;
                prd  = (gjit >= dm_i);
            else
                prd  = rand(numel(rb),1) < pw(rb);
            end

            mB = compute_metrics_soft(prd, y_b);
            BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
            BR.G(b)   = sqrt(max(mB.REC,0)*max(mB.SPEC,0));

            % regional bootstrap stats (soft)
            mask_ISPB_b = mask_ISPB_V(rb);
            mask_OSPB_b = mask_OSPB_V(rb);

            denom_w = sum(y_b(mask_ISPB_b), 'omitnan');
            denom_d = sum((1-y_b(mask_OSPB_b)), 'omitnan');

            if denom_w > 0
                BR.REC_ISPB_wet(b) = sum(double(prd(mask_ISPB_b)).*y_b(mask_ISPB_b), 'omitnan') / denom_w;
            end
            if denom_d > 0
                BR.REC_OSPB_dry(b) = sum(double(~prd(mask_OSPB_b)).*(1-y_b(mask_OSPB_b)), 'omitnan') / denom_d;
            end
            if isfinite(BR.REC_ISPB_wet(b)) && isfinite(BR.REC_OSPB_dry(b))
                BR.G_ISPBwet_OSPBdry(b) = sqrt(max(BR.REC_ISPB_wet(b),0) * max(BR.REC_OSPB_dry(b),0));
            end
        end

    else
        % Component bootstrap
        Vmask = false(numel(pw_all),1); Vmask(V) = true;
        GV = {};
        if ~isempty(G_M)
            GV = cellfun(@(g) g(Vmask(g)), G_M, 'uni', 0);
            GV = GV(~cellfun(@isempty, GV));
        end
        nC = numel(GV);

        if nC==0
            warning('[bootstrap] No valid components in ROI -> pixel bootstrap fallback');
            for b = 1:B
                rb = randi(Nf, Nf, 1);
                y_b = yr(rb);
                prd = rand(numel(rb),1) < pw(rb);

                mB = compute_metrics_soft(prd, y_b);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0)*max(mB.SPEC,0));

                mask_ISPB_b = mask_ISPB_V(rb);
                mask_OSPB_b = mask_OSPB_V(rb);

                denom_w = sum(y_b(mask_ISPB_b), 'omitnan');
                denom_d = sum((1-y_b(mask_OSPB_b)), 'omitnan');

                if denom_w > 0
                    BR.REC_ISPB_wet(b) = sum(double(prd(mask_ISPB_b)).*y_b(mask_ISPB_b), 'omitnan') / denom_w;
                end
                if denom_d > 0
                    BR.REC_OSPB_dry(b) = sum(double(~prd(mask_OSPB_b)).*(1-y_b(mask_OSPB_b)), 'omitnan') / denom_d;
                end
                if isfinite(BR.REC_ISPB_wet(b)) && isfinite(BR.REC_OSPB_dry(b))
                    BR.G_ISPBwet_OSPBdry(b) = sqrt(max(BR.REC_ISPB_wet(b),0) * max(BR.REC_OSPB_dry(b),0));
                end
            end
        else
            mfrac = max(0, min(1, cfg.uncertainty.mfrac));
            mapM2V = zeros(numel(Vmask),1); mapM2V(V) = 1:nnz(V);

            if cfg.uncertainty.pps_components
                wC = cellfun(@numel, GV); sW = sum(wC);
                if sW>0, wC = wC / sW; else, wC = []; end
            else
                wC = [];
            end

            for b = 1:B
                mC = max(1, round(mfrac * nC));
                if isempty(wC)
                    rbC = randi(nC, mC, 1);
                else
                    u = rand(mC,1); cdf = cumsum(wC(:));
                    rbC = arrayfun(@(x) find(cdf>=x,1,'first'), u);
                end

                pick = vertcat(GV{rbC});
                idxV = mapM2V(pick);

                pv = pw(idxV);
                yv = yr(idxV);

                prd = rand(numel(pv),1) < pv;

                mB = compute_metrics_soft(prd, yv);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0)*max(mB.SPEC,0));

                mask_ISPB_b = mask_ISPB_V(idxV);
                mask_OSPB_b = mask_OSPB_V(idxV);

                denom_w = sum(yv(mask_ISPB_b), 'omitnan');
                denom_d = sum((1-yv(mask_OSPB_b)), 'omitnan');

                if denom_w > 0
                    BR.REC_ISPB_wet(b) = sum(double(prd(mask_ISPB_b)).*yv(mask_ISPB_b), 'omitnan') / denom_w;
                end
                if denom_d > 0
                    BR.REC_OSPB_dry(b) = sum(double(~prd(mask_OSPB_b)).*(1-yv(mask_OSPB_b)), 'omitnan') / denom_d;
                end
                if isfinite(BR.REC_ISPB_wet(b)) && isfinite(BR.REC_OSPB_dry(b))
                    BR.G_ISPBwet_OSPBdry(b) = sqrt(max(BR.REC_ISPB_wet(b),0) * max(BR.REC_OSPB_dry(b),0));
                end
            end
        end
    end

    % CIs
    for s = stat_list
        vR = BR.(s{1});
        CI_raw.(s{1})(i,:)  = prctile(vR,  100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
        CIi_raw.(s{1})(i,:) = prctile(vR,  100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    end
    CI_raw.G(i,:)  = prctile(BR.G,  100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CIi_raw.G(i,:) = prctile(BR.G,  100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);

    CI_raw.REC_ISPB_wet(i,:)      = prctile(BR.REC_ISPB_wet,      100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CI_raw.REC_OSPB_dry(i,:)      = prctile(BR.REC_OSPB_dry,      100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CI_raw.G_ISPBwet_OSPBdry(i,:) = prctile(BR.G_ISPBwet_OSPBdry, 100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);

    CIi_raw.REC_ISPB_wet(i,:)      = prctile(BR.REC_ISPB_wet,      100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    CIi_raw.REC_OSPB_dry(i,:)      = prctile(BR.REC_OSPB_dry,      100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    CIi_raw.G_ISPBwet_OSPBdry(i,:) = prctile(BR.G_ISPBwet_OSPBdry, 100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
end

% Suggested margins (same idea as original)
is_flat = startsWith(S.names,'Flat_')'; nf = find(~is_flat);
if ~isempty(nf)
    S.suggested_margin_mean   = median(S.sigma_gdif_mean(nf),   'omitnan');
    S.suggested_margin_median = median(S.sigma_gdif_median(nf), 'omitnan');
else
    S.suggested_margin_mean   = median(S.sigma_gdif_mean,   'omitnan');
    S.suggested_margin_median = median(S.sigma_gdif_median, 'omitnan');
end
fprintf('[Margin] suggested decision margin ≈ mean: %.3f  | median: %.3f (mW m^-2)\n', ...
    S.suggested_margin_mean, S.suggested_margin_median);
fprintf('[Margin] mode: "%s" | fixed=%.3f\n', string(cfg.decision_margin_mode), cfg.decision_margin);

% Build Gstats in the SAME field names as original (so plotting code works)
S.Gstats = struct();
for s = stat_list
    S.Gstats.([s{1} '_raw']) = struct('EXP', EXP_raw.(s{1}), 'CI', CI_raw.(s{1}), 'CIi', CIi_raw.(s{1}));
end
S.Gstats.G_raw = struct('EXP', EXP_raw.G, 'CI', CI_raw.G, 'CIi', CIi_raw.G);

S.Gstats.REC_ISPB_wet_raw = struct('EXP', EXP_raw.REC_ISPB_wet, 'CI', CI_raw.REC_ISPB_wet, 'CIi', CIi_raw.REC_ISPB_wet);
S.Gstats.REC_OSPB_dry_raw = struct('EXP', EXP_raw.REC_OSPB_dry, 'CI', CI_raw.REC_OSPB_dry, 'CIi', CIi_raw.REC_OSPB_dry);
S.Gstats.G_ISPBwet_OSPBdry_raw = struct('EXP', EXP_raw.G_ISPBwet_OSPBdry, 'CI', CI_raw.G_ISPBwet_OSPBdry, 'CIi', CIi_raw.G_ISPBwet_OSPBdry);

S.Gstats.truth_mode = S.truth_mode_used;

%% ==================== Write CI summary CSV (like original) ====================
okCI = isfield(S,'Gstats') && isfield(S.Gstats,'SPEC_raw') && isfield(S.Gstats.SPEC_raw,'CI') && ...
       ~isempty(S.Gstats.SPEC_raw.CI) && size(S.Gstats.SPEC_raw.CI,2)==2;

if okCI
  try
    Name = string(S.names(:));
    if isfield(S,'titles') && numel(S.titles)==numel(S.names)
        Title = string(S.titles(:));
    else
        Title = Name;
    end

    n = numel(Name);
    getCol = @(fld, sub) i_getCol(S.Gstats, fld, sub, n);
    getCI  = @(fld, which) i_getCI(S.Gstats, fld, which, n);

    SPEC    = getCol('SPEC_raw','EXP');   SPEC_lo = getCI('SPEC_raw',1); SPEC_hi = getCI('SPEC_raw',2);
    REC     = getCol('REC_raw','EXP');    REC_lo  = getCI('REC_raw',1);  REC_hi  = getCI('REC_raw',2);
    PREC    = getCol('PREC_raw','EXP');   PREC_lo = getCI('PREC_raw',1); PREC_hi = getCI('PREC_raw',2);
    ACC     = getCol('ACC_raw','EXP');    ACC_lo  = getCI('ACC_raw',1);  ACC_hi  = getCI('ACC_raw',2);
    F1      = getCol('F1_raw','EXP');     F1_lo   = getCI('F1_raw',1);   F1_hi   = getCI('F1_raw',2);

    Tsum = table(Name, Title, ...
        SPEC, SPEC_lo, SPEC_hi, ...
        REC,  REC_lo,  REC_hi, ...
        PREC, PREC_lo, PREC_hi, ...
        ACC,  ACC_lo,  ACC_hi, ...
        F1,   F1_lo,   F1_hi, ...
        'VariableNames', {'Name','Title', ...
          'SPEC_raw','SPEC_raw_lo','SPEC_raw_hi', ...
          'REC_raw','REC_raw_lo','REC_raw_hi', ...
          'PREC_raw','PREC_raw_lo','PREC_raw_hi', ...
          'ACC_raw','ACC_raw_lo','ACC_raw_hi', ...
          'F1_raw','F1_raw_lo','F1_raw_hi'});

    if isfield(S.Gstats,'G_raw')
      G    = i_pad(S.Gstats.G_raw.EXP(:),  n);
      G_lo = i_pad(S.Gstats.G_raw.CI(:,1), n);
      G_hi = i_pad(S.Gstats.G_raw.CI(:,2), n);
      Tsum.G_raw    = G;
      Tsum.G_raw_lo = G_lo;
      Tsum.G_raw_hi = G_hi;
    end

    % Regional columns if present
    if isfield(S.Gstats,'REC_ISPB_wet_raw') && isfield(S.Gstats,'REC_OSPB_dry_raw') && isfield(S.Gstats,'G_ISPBwet_OSPBdry_raw')
        Tsum.REC_ISPB_wet_raw    = i_pad(S.Gstats.REC_ISPB_wet_raw.EXP(:), n);
        Tsum.REC_ISPB_wet_raw_lo = i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,1), n);
        Tsum.REC_ISPB_wet_raw_hi = i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,2), n);

        Tsum.REC_OSPB_dry_raw    = i_pad(S.Gstats.REC_OSPB_dry_raw.EXP(:), n);
        Tsum.REC_OSPB_dry_raw_lo = i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,1), n);
        Tsum.REC_OSPB_dry_raw_hi = i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,2), n);

        Tsum.G_ISPBwet_OSPBdry_raw    = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.EXP(:), n);
        Tsum.G_ISPBwet_OSPBdry_raw_lo = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,1), n);
        Tsum.G_ISPBwet_OSPBdry_raw_hi = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,2), n);
    end

    timestamp_ci = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    csv_path = fullfile(cfg.outdir, sprintf('ci_summary_%s_%s_%s.csv', region_suffix(cfg), S.truth_mode_used, timestamp_ci));
    writetable(Tsum, csv_path);
    fprintf('Saved CI summary CSV: %s\n', csv_path);

  catch
    warning('CI/diagnostic readout failed');
  end
end

%% ==================== Save lightweight plotting bundle (like original) =====
cfg.artifact_dir    = fullfile(cfg.outdir, 'artifacts');
cfg.ds_plot         = 8;
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end

ds = max(1, round(cfg.ds_plot));
r = 1:ds:size(S.H,1); c = 1:ds:size(S.H,2);

S_plot = struct();
S_plot.includes = {'X_km','Y_km','H','viz_mask','sink_mask_comp','comp_id', ...
                'roi_mask','y_raw_full','y_bin_full','y_w_full','truth_mode_used', ...
                'mask_ISPB','mask_OSPB','mask_SPB','rect_mask','ds_stride','names','titles','Gstats'};

% Downsampled model grids (optional)
S_plot.models_ds = struct();
for ii = 1:numel(S.names)
    fld = validNames{ii};
    if isfield(S.models, fld)
        S_plot.models_ds.(fld) = single(S.models.(fld)(r,c));
    end
end

S_plot.names  = S.names;
if isfield(S,'titles') && numel(S.titles)==numel(S.names), S_plot.titles = S.titles; else, S_plot.titles = S.names; end
S_plot.run_id = cfg.run_id;

S_plot.Gstats = S.Gstats;
S_plot.truth_mode_used = string(S.truth_mode_used);
S_plot.model_cvals = S.model_cvals;

S_plot.sigma_gdif_mean   = S.sigma_gdif_mean;
S_plot.sigma_gdif_median = S.sigma_gdif_median;
S_plot.suggested_margin_mean   = S.suggested_margin_mean;
S_plot.suggested_margin_median = S.suggested_margin_median;

S_plot.X_km     = single(S.X_km(r,c));
S_plot.Y_km     = single(S.Y_km(r,c));
S_plot.H        = single(S.H(r,c));
S_plot.viz_mask = logical(S.viz_mask(r,c));

S_plot.sink_mask_comp = logical(S.sink_mask_comp(r,c));
S_plot.comp_id        = uint32(S.comp_id(r,c));
S_plot.roi_mask       = logical(S.roi_mask(r,c));

S_plot.y_raw_full = single(S.y_raw_full(r,c));
S_plot.y_bin_full = single(S.y_bin_full(r,c));
S_plot.y_w_full   = single(S.y_w_full(r,c));

if isfield(S,'mask_ISPB'), S_plot.mask_ISPB = logical(S.mask_ISPB(r,c)); end
if isfield(S,'mask_OSPB'), S_plot.mask_OSPB = logical(S.mask_OSPB(r,c)); end
if isfield(S,'mask_SPB'),  S_plot.mask_SPB  = logical(S.mask_SPB(r,c));  end
if isfield(S,'rect_mask'), S_plot.rect_mask = logical(S.rect_mask(r,c)); end

S_plot.ds_stride = uint16(ds);
S_plot.margin_used_mode = string(cfg.decision_margin_mode);

cfg_plot = struct();
cfg_plot.outdir    = cfg.outdir;
cfg_plot.run_id    = cfg.run_id;
cfg_plot.font_size = cfg.font_size;
cfg_plot.overlay_alpha = cfg.overlay_alpha;
cfg_plot.truth = cfg_used.truth;
if isfield(cfg,'region'), cfg_plot.region = cfg.region; else, cfg_plot.region = struct('mode','ALL'); end
cfg_plot.flat = struct('interp_enable',true,'interp_points',50,'interp_method','spline');

stamp  = char(datetime('now','Format','yyyyMMdd_HHmmss'));
f_plot = fullfile(cfg.artifact_dir, sprintf('S_plot_small_%s_%s.mat', S.truth_mode_used, stamp));
save(f_plot, 'S_plot', 'cfg_plot', '-v7');
fprintf('[artifacts] Saved lightweight plot file: %s\n', f_plot);

%% ==================== Cache full-grid GΔ for non-flat models (like original) ==
cache_dir = fullfile(cfg.outdir, 'gdif_cache');
if ~exist(cache_dir,'dir'), mkdir(cache_dir); end

is_flat = startsWith(S.names,'Flat_');
for ii = 1:numel(S.names)
    if is_flat(ii), continue; end
    name = S.names{ii};
    fld  = validNames{ii};
    zd_path = fullfile(cache_dir, sprintf('Gdif_full_%s.mat', name));
    if exist(zd_path,'file'), continue; end

    if isfield(S,'Gdiff') && isfield(S.Gdiff,fld)
        Gdif_full = single(S.Gdiff.(fld));
    else
        [~, ~, tmpG] = processGHF(S.models.(fld), S.Ts, S.Mb, S.H);
        if isfield(S,'Delta_q_adv') && ~isempty(S.Delta_q_adv), tmpG = tmpG + S.Delta_q_adv; end
        Gdif_full = single(tmpG);
    end

    try
        save(zd_path, 'Gdif_full', '-v7');
        fprintf('[artifacts] cached GΔ -> %s\n', zd_path);
    catch ME
        warning('Failed to cache GΔ for %s: %s', name, ME.message);
    end
end

fprintf('\nDone. S contains:\n');
fprintf('  - spec_stats (ROI/VIZ/REG)\n');
fprintf('  - y_bin_vec / y_w_vec / y_use_vec (and *_full grids)\n');
fprintf('  - results_table_binary (legacy)\n');
fprintf('  - results_table_weighted (soft point metrics)\n');
fprintf('  - Gstats (truth_mode_used), CI summary CSV, S_plot_small, gdif_cache\n');

%% ==================== WEIGHTING IMPACT DIAGNOSTICS (PASTE BLOCK) ====================
% Purpose:
%   - Quantify + visualize how weighted truth (y_w) differs from binary truth (y_bin)
%   - Show how model ranking/score changes when using weighted vs binary point metrics
% Primary model score (per your definition):
%   G = REC * SPEC  (ROI-wide wet recall times dry recall)
%
% Assumes these already exist ABOVE this block:
%   Sc_roi, S.y_bin_vec, S.y_w_vec, cfg, cfg_used, S.truth_mode_used, S.names
%   S.results_table_binary, S.results_table_weighted
%   S.X_km, S.Y_km, S.y_bin_full, S.y_w_full, cfg.run_id, cfg.outdir, cfg.region.mode
%   (Optional) S.truth_cal

try
    diagdir = fullfile(cfg.outdir, 'weight_diagnostics');
    if ~exist(diagdir,'dir'), mkdir(diagdir); end

    %% ---- Basic ROI label stats ----
    Scv = double(Sc_roi(:));
    Scv = Scv(isfinite(Scv) & Scv~=9999);

    yB = double(S.y_bin_vec(:));
    yW = double(S.y_w_vec(:));
    yW = yW(isfinite(yW)); % should already be finite where ROI is valid

    fprintf('\n====== WEIGHT DIAGNOSTICS (ROI) ======\n');
    fprintf('region=%s | truth_mode_used=%s\n', string(cfg.region.mode), string(S.truth_mode_used));
    fprintf('ROI N=%d | Sc valid N=%d\n', numel(S.y_bin_vec), numel(Scv));
    if ~isempty(Scv)
        fprintf('Sc quantiles: p05=%.3f p25=%.3f p50=%.3f p75=%.3f p95=%.3f\n', ...
            quantile(Scv,0.05), quantile(Scv,0.25), quantile(Scv,0.50), quantile(Scv,0.75), quantile(Scv,0.95));
    end
    fprintf('binary wet fraction (Sc>%.2f): %.3f\n', cfg.spec_thresh, mean(yB>0));
    fprintf('weighted mean(E[y])=%.3f | median=%.3f\n', mean(yW,'omitnan'), median(yW,'omitnan'));

    bins = [0 0.1 0.25 0.5 0.75 0.9 1.0];
    h = histcounts(yW, bins, 'Normalization','probability');
    fprintf('weight mass bins [0-0.1,0.1-0.25,0.25-0.5,0.5-0.75,0.75-0.9,0.9-1]:\n');
    fprintf('  %s\n', sprintf('%.3f ', h));

    fprintf('disagreement(binary vs weight>0.5): %.3f\n', mean((yB>0) ~= (yW>0.5)));

    % Approx effective sample size (ESS) for wet/dry mass (informal but intuitive)
    w  = yW(:); w  = w(isfinite(w));
    w0 = 1 - w;
    ESS_wet = (sum(w)^2)  / max(sum(w.^2),  eps);
    ESS_dry = (sum(w0)^2) / max(sum(w0.^2), eps);
    fprintf('approx ESS wet=%.1f | ESS dry=%.1f (out of N=%d)\n', ESS_wet, ESS_dry, numel(yW));

    if isfield(S,'truth_cal')
        tc = S.truth_cal;
        fprintf('logistic auto-cal anchors: Sc_low(q=%.2f)=%.3f -> w=%.2f | Sc_high(q=%.2f)=%.3f -> w=%.2f\n', ...
            tc.q_low, tc.Sc_low, tc.w_low, tc.q_high, tc.Sc_high, tc.w_high);
        fprintf('logistic params: tmid=%.3f  k=%.3f\n', tc.tmid, tc.k);
    end
    fprintf('=====================================\n\n');

    %% ---- Figure 1: Sc histogram + threshold ----
    f1 = figure('Color','w');
    histogram(Scv, 60, 'Normalization','pdf'); hold on;
    xline(cfg.spec_thresh, '--', sprintf('thresh=%.2f', cfg.spec_thresh));
    xlabel('Specularity content (Sc)'); ylabel('PDF');
    title(sprintf('ROI Sc distribution (%s)', string(cfg.region.mode)));
    grid on;
    saveas(f1, fullfile(diagdir, sprintf('Sc_hist_%s.png', cfg.run_id)));

    %% ---- Figure 2: weight function w(Sc) + anchors ----
    if ~isempty(Scv)
        sgrid = linspace(min(Scv), max(Scv), 300);
    else
        sgrid = linspace(0, 1, 300);
    end
    wgrid = spec_weight(sgrid, cfg_used.truth);

    f2 = figure('Color','w');
    plot(sgrid, wgrid, 'LineWidth', 2); hold on;
    xline(cfg.spec_thresh, '--', 'binary thresh');
    yline(0.5, ':');
    xlabel('Sc'); ylabel('weight w(Sc)');
    title('Weighting function');
    grid on;
    if isfield(S,'truth_cal')
        tc = S.truth_cal;
        plot([tc.Sc_low tc.Sc_high], [tc.w_low tc.w_high], 'o', 'MarkerSize',8, 'LineWidth',2);
        xline(tc.tmid, ':', sprintf('tmid=%.3f', tc.tmid));
    end
    saveas(f2, fullfile(diagdir, sprintf('weight_function_%s.png', cfg.run_id)));

    %% ---- Figure 3: ROI scatter Sc vs weight ----
    % (use Sc_roi directly to preserve pairing)
    vmask = isfinite(Sc_roi) & (Sc_roi~=9999) & isfinite(S.y_w_vec);
    f3 = figure('Color','w');
    scatter(double(Sc_roi(vmask)), double(S.y_w_vec(vmask)), 8, 'filled'); hold on;
    xline(cfg.spec_thresh, '--');
    yline(0.5, ':');
    xlabel('Sc (ROI)'); ylabel('weight');
    title('ROI: Sc vs weight');
    grid on;
    saveas(f3, fullfile(diagdir, sprintf('Sc_vs_weight_%s.png', cfg.run_id)));

    %% ---- Figure 4: maps of weight + (weight - binary) + ambiguous zone ----
    % Map: weight
    if isfield(S,'y_w_full') && isfield(S,'y_bin_full') && isfield(S,'X_km') && isfield(S,'Y_km')
        f4 = figure('Color','w');
        imagesc(S.X_km(1,:), S.Y_km(:,1), S.y_w_full); axis image; set(gca,'YDir','normal');
        title(sprintf('Weighted truth w (full grid; ROI non-NaN) | %s', string(cfg.region.mode)));
        xlabel('X (km)'); ylabel('Y (km)'); colorbar;
        saveas(f4, fullfile(diagdir, sprintf('map_weight_%s.png', cfg.run_id)));

        % Map: difference from binary
        D = double(S.y_w_full) - double(S.y_bin_full);
        f5 = figure('Color','w');
        imagesc(S.X_km(1,:), S.Y_km(:,1), D); axis image; set(gca,'YDir','normal');
        title('Weighting effect map: w - y_{bin}');
        xlabel('X (km)'); ylabel('Y (km)'); colorbar;
        saveas(f5, fullfile(diagdir, sprintf('map_weight_minus_binary_%s.png', cfg.run_id)));

        % Map: ambiguous zone
        amb = (S.y_w_full >= 0.25) & (S.y_w_full <= 0.75);
        f6 = figure('Color','w');
        imagesc(S.X_km(1,:), S.Y_km(:,1), amb); axis image; set(gca,'YDir','normal');
        title('Ambiguous-weight zone (0.25 \le w \le 0.75)');
        xlabel('X (km)'); ylabel('Y (km)'); colorbar;
        saveas(f6, fullfile(diagdir, sprintf('map_ambiguous_zone_%s.png', cfg.run_id)));
    end

    %% ---- Model-score impact (your primary score: G = REC * SPEC) ----
    Tb = S.results_table_binary;
    Tw = S.results_table_weighted;

    % Binary ROI score
    Gb = double(Tb.REC) .* double(Tb.SPEC);

    % Weighted ROI score (soft metrics)
    % (your Tw uses RECw and SPECw as the ROI-wide wet/dry recalls under soft truth)
    Gw = double(Tw.RECw) .* double(Tw.SPECw);

    % Rank (1=best)
    [~,ordb] = sort(Gb,'descend');
    [~,ordw] = sort(Gw,'descend');
    rank_b = zeros(numel(Gb),1); rank_b(ordb) = 1:numel(Gb);
    rank_w = zeros(numel(Gw),1); rank_w(ordw) = 1:numel(Gw);

    dG    = Gw - Gb;
    dRank = rank_w - rank_b;

    Tcmp = table(string(Tb.Model), Gb, Gw, dG, rank_b, rank_w, dRank, ...
        'VariableNames', {'Model','G_binary','G_weighted','dG','rank_binary','rank_weighted','dRank'});

    cmp_csv = fullfile(diagdir, sprintf('G_rank_compare_%s.csv', cfg.run_id));
    writetable(Tcmp, cmp_csv);
    fprintf('[weight diagnostics] wrote: %s\n', cmp_csv);

    % Scatter: G_binary vs G_weighted (quick visual of changes)
    f7 = figure('Color','w');
    scatter(Gb, Gw, 40, 'filled'); grid on;
    xlabel('G = REC*SPEC (binary point metrics)'); ylabel('G = RECw*SPECw (weighted point metrics)');
    title('Model score comparison: binary vs weighted');
    saveas(f7, fullfile(diagdir, sprintf('G_scatter_binary_vs_weighted_%s.png', cfg.run_id)));

    % Bar plot: rank shift
    f8 = figure('Color','w');
    bar(dRank); grid on;
    xlabel('Model index (as in S.names)'); ylabel('\Delta rank (weighted - binary)');
    title('\Delta rank from weighting (negative = improves under weighting)');
    saveas(f8, fullfile(diagdir, sprintf('rank_shift_%s.png', cfg.run_id)));

catch ME
    fprintf('Weight diagnostics failed');
end

end % ===== function end =====

%% ===================== Helpers =====================

function st = summarize_SQ(Q, mask)
    q = Q(mask);
    q = q(isfinite(q) & q~=9999);
    st = struct('N',0,'min',NaN,'p01',NaN,'p05',NaN,'p10',NaN,'p25',NaN,'p50',NaN,'p75',NaN,'p90',NaN,'p95',NaN,'p99',NaN,'max',NaN,'mean',NaN,'std',NaN);
    if isempty(q), return; end
    st.N    = numel(q);
    st.min  = min(q);
    st.max  = max(q);
    st.mean = mean(q);
    st.std  = std(q);
    p = [0.01 0.05 0.10 0.25 0.50 0.75 0.90 0.95 0.99];
    v = quantile(q, p);
    st.p01=v(1); st.p05=v(2); st.p10=v(3); st.p25=v(4); st.p50=v(5);
    st.p75=v(6); st.p90=v(7); st.p95=v(8); st.p99=v(9);
end

function w = spec_weight(Sc, truthCfg)
    w = zeros(size(Sc), 'double');
    valid = isfinite(Sc) & (Sc~=9999);
    if ~any(valid), return; end

    mode = lower(string(truthCfg.weight_mode));
    s = double(Sc(valid));

    switch mode
        case "ramp"
            t0 = truthCfg.t0;
            wv = (s - t0) ./ max(1 - t0, eps);
            wv = min(1, max(0, wv));

        case "logistic"
            tmid = truthCfg.tmid;
            k    = truthCfg.k;
            wv = 1 ./ (1 + exp(-k*(s - tmid)));

        case "power"
            t0 = truthCfg.t0;
            p  = truthCfg.p;
            wv = (s - t0) ./ max(1 - t0, eps);
            wv = min(1, max(0, wv));
            wv = wv.^p;

        otherwise
            error('Unknown truth.weight_mode: %s', truthCfg.weight_mode);
    end

    w(valid) = wv;
end

function m = compute_metrics_hard(pred, truth)
    pred  = logical(pred(:));
    truth = logical(truth(:));
    N  = numel(truth);
    TP = sum( pred &  truth);
    FP = sum( pred & ~truth);
    TN = sum(~pred & ~truth);
    FN = sum(~pred &  truth);
    ACC  = (TP + TN) / max(N, eps);
    PREC = TP / max(TP + FP, eps);
    SPEC = TN / max(TN + FP, eps);
    REC  = TP / max(TP + FN, eps);
    F1   = 2 * PREC * REC / max(PREC + REC, eps);
    m = struct('N',N,'TP',TP,'FP',FP,'TN',TN,'FN',FN,'ACC',ACC,'PREC',PREC,'SPEC',SPEC,'REC',REC,'F1',F1);
end

function m = compute_metrics_soft(pred, truth_w)
    pred = double(logical(pred(:)));
    y    = double(truth_w(:));
    y    = min(1, max(0, y));

    TP = sum(pred .* y, 'omitnan');
    FP = sum(pred .* (1 - y), 'omitnan');
    FN = sum((1 - pred) .* y, 'omitnan');
    TN = sum((1 - pred) .* (1 - y), 'omitnan');

    N = TP + FP + FN + TN;
    ACC  = (TP + TN) / max(N, eps);
    PREC = TP / max(TP + FP, eps);
    SPEC = TN / max(TN + FP, eps);
    REC  = TP / max(TP + FN, eps);
    F1   = 2 * PREC * REC / max(PREC + REC, eps);

    m = struct('N',N,'TP',TP,'FP',FP,'TN',TN,'FN',FN,'ACC',ACC,'PREC',PREC,'SPEC',SPEC,'REC',REC,'F1',F1);
end

function pwet = wet_prob_with_margin(gdif_vec, sig_m, sig_g, margin)
  if isempty(sig_g), sig_g = 0; end
  if nargin < 4 || isempty(margin), margin = 0; end
  s2 = sig_m.^2 + sig_g.^2;
  pwet = NaN(size(gdif_vec));
  Z = (s2==0);
  pwet(Z) = double(gdif_vec(Z) >= margin);
  NZ = ~Z;
  if any(NZ(:))
      s = sqrt(max(s2(NZ), eps));
      pwet(NZ) = 0.5*(1 + erf(((gdif_vec(NZ) - margin)./s)/sqrt(2)));
  end
end

function out = merge_structs(base, override)
    out = base;
    fn = fieldnames(override);
    for i = 1:numel(fn)
        if isstruct(override.(fn{i})) && isfield(out, fn{i}) && isstruct(out.(fn{i}))
            out.(fn{i}) = merge_structs(out.(fn{i}), override.(fn{i}));
        else
            out.(fn{i}) = override.(fn{i});
        end
    end
end

function dm = get_model_margin(cfg, S, iModel)
    mode = 'fixed';
    if isfield(cfg,'decision_margin_mode') && ~isempty(cfg.decision_margin_mode)
        mode = lower(string(cfg.decision_margin_mode));
    end

    haveMean   = isfield(S,'sigma_gdif_mean')   && ~isempty(S.sigma_gdif_mean);
    haveMedian = isfield(S,'sigma_gdif_median') && ~isempty(S.sigma_gdif_median);

    switch mode
        case "auto_mean"
            if isfield(S,'suggested_margin_mean') && isfinite(S.suggested_margin_mean)
                dm = S.suggested_margin_mean;
            elseif haveMean
                dm = median(S.sigma_gdif_mean, 'omitnan');
            else
                dm = cfg.decision_margin;
            end
        case "auto_median"
            if isfield(S,'suggested_margin_median') && isfinite(S.suggested_margin_median)
                dm = S.suggested_margin_median;
            elseif haveMedian
                dm = median(S.sigma_gdif_median, 'omitnan');
            else
                dm = cfg.decision_margin;
            end
        case "auto_per_model_mean"
            if haveMean && iModel>=1 && iModel<=numel(S.sigma_gdif_mean) && isfinite(S.sigma_gdif_mean(iModel))
                dm = S.sigma_gdif_mean(iModel);
            else
                dm = cfg.decision_margin;
            end
        case "auto_per_model_median"
            if haveMedian && iModel>=1 && iModel<=numel(S.sigma_gdif_median) && isfinite(S.sigma_gdif_median(iModel))
                dm = S.sigma_gdif_median(iModel);
            else
                dm = cfg.decision_margin;
            end
        otherwise
            dm = cfg.decision_margin;
    end
    if ~isfinite(dm) || isempty(dm), dm = cfg.decision_margin; end
end

function v = getfield_default(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end

function s = region_suffix(cfg)
    s = '';
    try
        if isfield(cfg,'region') && isfield(cfg.region,'mode') && ~isempty(cfg.region.mode)
            s = ['_' upper(char(cfg.region.mode))];
        end
    catch
    end
end

function v = i_pad(v, n)
  v = v(:);
  if numel(v) < n, v(end+1:n,1) = NaN; end
  if numel(v) > n, v = v(1:n); end
end
function v = i_getCol(GS, field, sub, n)
  v = NaN(n,1);
  if isfield(GS,field) && isfield(GS.(field),sub) && ~isempty(GS.(field).(sub))
    v = i_pad(GS.(field).(sub)(:), n);
  end
end
function v = i_getCI(GS, field, which, n)
  v = NaN(n,1);
  if isfield(GS,field) && isfield(GS.(field),'CI') && ~isempty(GS.(field).CI) && size(GS.(field).CI,2)>=which
    v = i_pad(GS.(field).CI(:,which), n);
  end
end

function Mfull = expand_or_resize_mask(maskMaybe, meta, targetSize, name)
    nr = targetSize(1); nc = targetSize(2);
    Mfull = false(nr,nc);

    if isempty(maskMaybe)
        return
    end

    if isequal(size(maskMaybe), [nr nc])
        Mfull = logical(maskMaybe);
        return
    end

    haveCrop = all(isfield(meta, {'rmin','rmax','cmin','cmax'}));
    if haveCrop
        rr = meta.rmin:meta.rmax; cc = meta.cmin:meta.cmax;
        if numel(rr)==size(maskMaybe,1) && numel(cc)==size(maskMaybe,2)
            Mfull(rr,cc) = logical(maskMaybe);
            return
        end
    end

    warning('[mask] %s has size %dx%d; resizing to %dx%d (nearest).', ...
        name, size(maskMaybe,1), size(maskMaybe,2), nr, nc);
    Mfull = logical(imresize(maskMaybe, [nr nc], 'nearest'));
end
