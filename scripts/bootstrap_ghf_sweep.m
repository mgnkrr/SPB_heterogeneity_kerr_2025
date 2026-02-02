% ========================================================================
% bootstrap_ghf_sweep.m  (single-file, self-contained, functionality-preserving)
%
% PURPOSE
%   Run a sweep over specularity thresholds WITHOUT re-loading/re-building
%   static data each time.
%
% IMPORTANT CHANGE (per your request)
%   - REG_MASK IS REMOVED ENTIRELY.
%   - ROI/VIZ are now defined over the full slow+valid domain (optionally sink-gated).
%   - ISPB/OSPB masks are still loaded/stitched ONLY to support the regional metrics.
%
% USAGE
%   OUT = bootstrap_ghf_sweep(); % default cfg sweep
%   OUT = bootstrap_ghf_sweep(struct('sweep_mode',false,'spec_thresh',0.15));
%
% NOTE
%   Requires your project functions on path if datasets are missing:
%     build_data_core, add_g_fields, add_g_uncertainty, build_spb_masks, make_sink_mask
%
% ========================================================================

function OUT = bootstrap_ghf_sweep(cfg_override)

dbstop if error

% ---- ensure this script folder + datasets are on path ----
this_dir = fileparts(mfilename('fullpath'));
addpath(this_dir);
addpath(fullfile(this_dir,'datasets'));

% ---- deps once (best-effort) ----
check_dependencies_local(struct('autoInstall',true,'needAMT',false,'verbose',true));

% ---- cfg ----
cfg = make_cfg();
if nargin > 0 && ~isempty(cfg_override)
    cfg = merge_structs(cfg, cfg_override);
end

if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
cfg.artifact_dir = fullfile(cfg.outdir,'artifacts');
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end

% ---- prepare static once ----
STATIC = prepare_static(cfg);

% ---- run ----
if isfield(cfg,'sweep_mode') && cfg.sweep_mode
    OUT = run_sweep(cfg, STATIC);
else
    OUT = run_once(cfg, STATIC);
end

end

%% ========================================================================
%  CONFIG
% ========================================================================
function cfg = make_cfg()

run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));

cfg = struct( ...
  'spec_thresh',         0.2, ...    % Specularity content threshold
  'v_keep',              10,  ...    % Ice velocity threshold (m/yr)
  'decision_margin',     0,   ...    % Fixed margin: wet if (GΔ >= margin)
  'decision_margin_mode','fixed', ...% fixed|auto_mean|auto_median|auto_per_model_mean|auto_per_model_median
  'outdir',              'figs_out', ...
  'overwrite',           true, ...
  'overlay_alpha',       0.4, ...
  'font_size',           20, ...
  'to_single',           true, ...
  'sweep_mode',          true ...
);

cfg.run_id = run_id;
cfg.save_artifacts = true;

% Evaluation controls
cfg.eval = struct( ...
  'use_sinks_in_roi', false ...      % if true, ROI = slow&valid & sink_mask_comp
);

% Sinks (viz)
cfg.sinks = struct('marker','.', 'size',4, 'color',[0.05 0.05 0.05], 'alpha',0.85);

% Uncertainty
cfg.uncertainty = struct( ...
  'mode','analytic', ...
  'bootstrap_mode','component', ...  % pixel | component
  'n_boot',10, ...
  'seed',42, ...
  'mfrac',0.7, ...
  'band_main',0.95, ...
  'band_inner',0.50, ...
  'noise_enable',false, ...
  'noise_K',1, ...
  'pred_mode','bernoulli', ...       % bernoulli | noisy_threshold
  'class_balance',false, ...         % pixel bootstrap: stratified by class
  'pps_components',false ...         % component bootstrap: PPS by component size
);

cfg.plot_labels = struct('avoid_overlap',true,'max_iter',200,'step',0.01,'pad',0.005);

% Sweep list default
cfg.sweep = struct('spec_list', (10:5:30)/100);

% Files (centralize)
cfg.files = struct();
cfg.files.gmin_data = '/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/gmin_data.mat';
cfg.files.adv_maps  = 'datasets/longitudinal_advection_maps.mat';
cfg.files.spb_masks = fullfile(fileparts(mfilename('fullpath')), 'datasets', 'spb_masks.mat');
cfg.files.sink_new  = '/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_new.mat';
cfg.files.sink_cmp  = '/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat';

end

%% ========================================================================
%  STATIC PREP (RUN ONCE)
% ========================================================================
function STATIC = prepare_static(cfg)

fprintf('Preparing static data (one-time) ...\n');

STATIC = struct();
STATIC.cfg_base = cfg;  % for record

% -------------------- Load gmin_data.mat --------------------
fprintf('[static] Loading data ...\n');
datafile = cfg.files.gmin_data;

if exist(datafile,'file') == 2
    fprintf('[static] Loading %s ... ', datafile);
    L = load(datafile,'S');
    S = L.S; clear L;
    fprintf('ok.\n');
else
    warning('Data file not found: %s\nAttempting to build now...', datafile);
    require_on_path({'build_data_core','add_g_fields','add_g_uncertainty'});
    S = build_data_core('to_single',true,'save',true, ...
                    'outdir','datasets','outfile','gmin_data.mat');
    S = add_g_fields(fullfile('datasets','gmin_data.mat'), ...
                 'advection_file',cfg.files.adv_maps, ...
                 'advection_var','Delta_q_adv', 'advection_units','mWm2', ...
                 'save',true, 'outfile', fullfile('datasets','gmin_data.mat'));
    S = add_g_uncertainty('datasets/gmin_data.mat', ...
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

fprintf('[static] Grid: %dx%d | models: %d\n', size(S.Xgrid,1), size(S.Xgrid,2), numel(S.names));

% -------------------- Load/build SPB masks (ONLY for ISPB/OSPB regional metrics) --------------------
require_on_path({'build_spb_masks'}); % may be needed if file missing

maskfile = cfg.files.spb_masks;
if exist(maskfile,'file') == 2
    m = load(maskfile);
else
    fprintf('[static] spb_masks.mat missing -> running build_spb_masks()...\n');
    old = pwd;
    cd(fileparts(maskfile));
    build_spb_masks();
    cd(old);
    if exist(maskfile,'file') ~= 2
        error('build_spb_masks() did not create expected file: %s', maskfile);
    end
    m = load(maskfile);
end

S = stitch_spb_masks_to_grid(S, m);  % builds mask_ISPB/mask_OSPB/mask_SPB on S grid

% -------------------- Load/build sink masks --------------------
require_on_path({'make_sink_mask'}); % used if sink files missing

sink_new = cfg.files.sink_new;
sink_cmp = cfg.files.sink_cmp;

if exist(sink_new,'file')==2 && exist(sink_cmp,'file')==2
    t = load(sink_new, 'sink_mask'); S.sink_mask = t.sink_mask; clear t;
    t = load(sink_cmp, 'comp_id','sink_mask_comp');
    if isfield(t,'comp_id'), S.comp_id = t.comp_id; end
    S.sink_mask_comp = t.sink_mask_comp; clear t;
else
    fprintf('[static] sink masks missing -> running make_sink_mask()...\n');
    make_sink_mask();
    t = load(sink_new, 'sink_mask'); S.sink_mask = t.sink_mask; clear t;
    t = load(sink_cmp, 'comp_id','sink_mask_comp');
    if isfield(t,'comp_id'), S.comp_id = t.comp_id; end
    S.sink_mask_comp = t.sink_mask_comp; clear t;
end

% Force sink masks to S grid if needed (nearest, logical-safe)
[S.sink_mask,      did1] = force_mask_to_size(S.sink_mask,      size(S.Xgrid), 'sink_mask');
[S.sink_mask_comp, did2] = force_mask_to_size(S.sink_mask_comp, size(S.Xgrid), 'sink_mask_comp');
if did1 || did2
    warning('Sink masks were resized to match S grid. Ideally regenerate them on this grid.');
end

% Ensure logical
S.sink_mask      = logical(S.sink_mask);
S.sink_mask_comp = logical(S.sink_mask_comp);

% If comp_id missing/mismatched, rebuild from sink_mask_comp
if ~isfield(S,'comp_id') || isempty(S.comp_id) || ~isequal(size(S.comp_id), size(S.sink_mask_comp))
    S.comp_id = bwlabel(S.sink_mask_comp, 8);
end

% -------------------- Build ROI / VIZ masks (static, NOT dependent on spec_thresh) --------------------
valid_mask = isfinite(S.Q) & ~S.spec_invalid & isfinite(S.H) & isfinite(S.icevel);
slow_mask  = valid_mask & ~(S.icevel > cfg.v_keep);

% NO REG_MASK: evaluate everywhere slow+valid
S.viz_mask = slow_mask;

% ROI: default = slow+valid; optionally gate by sinks
if isfield(cfg,'eval') && isfield(cfg.eval,'use_sinks_in_roi') && cfg.eval.use_sinks_in_roi
    S.roi_mask = S.viz_mask & S.sink_mask_comp;
else
    S.roi_mask = S.viz_mask;
end
S.roi_mask = logical(S.roi_mask);

idx = find(S.roi_mask);
S.eval_mask_full = S.roi_mask;

% M-space selector (you can later gate further)
S.eval_mask = true(numel(idx),1);

% Median-domain index (static)
med_mask = S.viz_mask;
idx_med  = find(med_mask);

% -------------------- Component groups for component bootstrap (static) --------------------
% IMPORTANT: groups are stored as ROI-vector indices (1..numel(idx)), not M-space.
S.boot_groups = struct('nC',0,'groups',{{}},'uniqC',[]);
if strcmpi(cfg.uncertainty.bootstrap_mode,'component')
    comp_roi = S.comp_id(idx);     % comp ids for ROI pixels
    keep     = comp_roi > 0;       % sink components only
    if any(keep)
        comp_keep = comp_roi(keep);
        roi_pos   = find(keep);    % ROI-vector indices
        [uniqC,~,gidx] = unique(comp_keep);
        groups = accumarray(gidx, roi_pos, [], @(v){v});
        S.boot_groups.uniqC  = uniqC;
        S.boot_groups.groups = groups;
        S.boot_groups.nC     = numel(groups);
    end
end

% -------------------- Cache per-model center values (static; uses idx_med) --------------------
nM = numel(S.names);
validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);

if ~isfield(S,'model_cvals'),        S.model_cvals = nan(nM,1);        end
if ~isfield(S,'model_cvals_median'), S.model_cvals_median = nan(nM,1); end

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};
    Mi     = S.models.(fld);
    if startsWith(name_i,'Flat_')
        v = sscanf(name_i,'Flat_%f');
        S.model_cvals_median(i) = v;
        S.model_cvals(i)        = v;
    else
        S.model_cvals_median(i) = median(Mi(idx_med), 'omitnan');
        S.model_cvals(i)        = S.model_cvals_median(i);
    end
end
S.model_quantiles = struct('q50', S.model_cvals_median(:));

% -------------------- Cache GΔ vectors in ROI space (static, huge win) --------------------
gdif_cache = cell(nM,1);
for i = 1:nM
    fld = validNames{i};
    gi  = S.Gdiff.(fld);
    gdif_cache{i} = gi(idx); % ROI-vector length
end

% -------------------- Cache region membership vectors in ROI space (static) --------------------
mask_ISPB_vec = false(numel(idx),1);
mask_OSPB_vec = false(numel(idx),1);
if isfield(S,'mask_ISPB') && ~isempty(S.mask_ISPB), mask_ISPB_vec = S.mask_ISPB(idx); end
if isfield(S,'mask_OSPB') && ~isempty(S.mask_OSPB), mask_OSPB_vec = S.mask_OSPB(idx); end

% -------------------- Cache σ(Gmin) vector in M-space (static) --------------------
M  = S.eval_mask;
NM = nnz(M);

if isfield(S,'sigma_gmin') && ~isempty(S.sigma_gmin) && isequal(size(S.sigma_gmin), size(S.H))
    sigG_roi = S.sigma_gmin(idx);
    sigG_M   = sigG_roi(M);
else
    sigG_M   = zeros(NM,1);
end

% -------------------- Scalar σ(model) fallback per model (static) --------------------
scalar_sigma_per_model = nan(nM,1);
sigma_gdif_mean   = nan(nM,1);
sigma_gdif_median = nan(nM,1);

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
        sigM_roi = sigM_grid(idx); sigM_M = sigM_roi(M);
        sigma_gdel_vec = sqrt(sigM_M.^2 + sigG_M.^2);
        sigma_gdif_mean(i)   = mean(sigma_gdel_vec, 'omitnan');
        sigma_gdif_median(i) = median(sigma_gdel_vec, 'omitnan');
        scalar_sigma_per_model(i) = mean(sigM_M, 'omitnan');
    else
        guess = NaN;
        if isfield(S,'sigma_model_scalar') && numel(S.sigma_model_scalar)>=i
            guess = S.sigma_model_scalar(i);
        elseif isfinite(sigma_gdif_mean(i))
            guess = sigma_gdif_mean(i);
        end
        if ~isfinite(guess), guess = 0; end
        scalar_sigma_per_model(i) = guess;

        sigma_gdif_mean(i)   = mean(sqrt((guess.^2) + (sigG_M.^2)), 'omitnan');
        sigma_gdif_median(i) = median(sqrt((guess.^2) + (sigG_M.^2)), 'omitnan');
    end
end

S.sigma_gdif_mean   = sigma_gdif_mean;
S.sigma_gdif_median = sigma_gdif_median;

% Suggested margins (static)
is_flat = startsWith(S.names,'Flat_')'; nf = find(~is_flat);
if ~isempty(nf)
    S.suggested_margin_mean   = median(S.sigma_gdif_mean(nf),   'omitnan');
    S.suggested_margin_median = median(S.sigma_gdif_median(nf), 'omitnan');
else
    S.suggested_margin_mean   = median(S.sigma_gdif_mean,   'omitnan');
    S.suggested_margin_median = median(S.sigma_gdif_median, 'omitnan');
end

% -------------------- Stash back --------------------
STATIC.S = S;
STATIC.idx = idx;
STATIC.idx_med = idx_med;
STATIC.M = M;
STATIC.NM = NM;
STATIC.validNames = validNames;
STATIC.gdif_cache = gdif_cache;
STATIC.mask_ISPB_vec = mask_ISPB_vec;
STATIC.mask_OSPB_vec = mask_OSPB_vec;
STATIC.sigG_M = sigG_M;
STATIC.scalar_sigma_per_model = scalar_sigma_per_model;

% Coverage print (static)
fprintf('\n---------- STATIC SANITY & COVERAGE ----------\n');
[nr,nc] = size(S.Xgrid);
fprintf('[grid] size=%dx%d  (%.2f Mpx)\n', nr, nc, numel(S.Xgrid)/1e6);
fprintf('[masks] VIZ=%d, ROI=%d  (ROI/VIZ=%.1f%%)\n', ...
    nnz(S.viz_mask), nnz(S.roi_mask), 100*nnz(S.roi_mask)/max(nnz(S.viz_mask),1));
fprintf('[ROI] use_sinks_in_roi=%d\n', isfield(cfg,'eval') && isfield(cfg.eval,'use_sinks_in_roi') && cfg.eval.use_sinks_in_roi);
fprintf('[components] nC=%d\n', S.boot_groups.nC);
fprintf('---------------------------------------------\n\n');

end

%% ========================================================================
%  SWEEP DRIVER (NO RELOADING)
% ========================================================================
function OUT = run_sweep(cfg, STATIC)

spec_list = cfg.sweep.spec_list;
base_outdir = cfg.outdir;

OUT = struct();
OUT.run_id = cfg.run_id;
OUT.spec_thresh_list = spec_list;
OUT.runs = cell(numel(spec_list),1);

for k = 1:numel(spec_list)
    st = spec_list(k);

    cfg_k = cfg;
    cfg_k.sweep_mode = false;
    cfg_k.spec_thresh = st;

    cfg_k.outdir = fullfile(base_outdir, sprintf('spec_%0.2f', st));
    if ~exist(cfg_k.outdir,'dir'), mkdir(cfg_k.outdir); end

    cfg_k.artifact_dir = fullfile(cfg_k.outdir,'artifacts');
    if ~exist(cfg_k.artifact_dir,'dir'), mkdir(cfg_k.artifact_dir); end

    cfg_k.run_id = sprintf('%s_spec%0.2f', cfg.run_id, st);

    fprintf('\n==================== SPEC SWEEP %d/%d: spec_thresh = %.2f ====================\n', ...
        k, numel(spec_list), st);

    OUT.runs{k} = run_once(cfg_k, STATIC);
end

end

%% ========================================================================
%  ONE RUN (PER SPEC_THRESH)
% ========================================================================
function OUT = run_once(cfg, STATIC)

S   = STATIC.S;               % local copy-on-write struct
idx = STATIC.idx;             % ROI-vector index into full grid
M   = STATIC.M;               % M-space selector (logical length numel(idx))
NM  = STATIC.NM;

nM = numel(S.names);
validNames = STATIC.validNames;
gdif_cache = STATIC.gdif_cache;

% -------------------- Labels for this threshold --------------------
spec_roi = double(S.Q(idx));
spec_ok  = isfinite(spec_roi) & (spec_roi ~= 9999) & (spec_roi > cfg.spec_thresh);

y_raw_vec = logical(spec_ok(:));     % ROI-vector space
y_raw_M   = double(y_raw_vec(M));    % M-space (eval subset)

fprintf('[diag] spec_thresh=%.2f | frac(wet in ROI)=%.4f | N=%d\n', cfg.spec_thresh, mean(y_raw_vec), numel(y_raw_vec));

% Optional: full-grid label field
if cfg.save_artifacts
    y_raw_full = false(size(S.Q));
    y_raw_full(idx) = y_raw_vec;
else
    y_raw_full = [];
end

% Regional masks in ROI-vector space
wet_vec = y_raw_vec;
dry_vec = ~wet_vec;

mask_wet_ISPB_vec = STATIC.mask_ISPB_vec & wet_vec;
mask_dry_OSPB_vec = STATIC.mask_OSPB_vec & dry_vec;

N_wet_ISPB = nnz(mask_wet_ISPB_vec);
N_dry_OSPB = nnz(mask_dry_OSPB_vec);

% Regional masks in M-space
mask_wet_ISPB_M = mask_wet_ISPB_vec(M);
mask_dry_OSPB_M = mask_dry_OSPB_vec(M);

% -------------------- Point metrics (no uncertainty) --------------------
fprintf('Evaluating point metrics on ROI...\n');

TP  = zeros(nM,1,'double'); FP  = TP; TN = TP; FN = TP;
ACC = TP; PR_r = TP; RC = TP; F1 = TP; TPFP=TP; TNFN=TP;
margin_used = nan(nM,1);

RAW_G = nan(nM,1);

REC_ISPB_wet      = nan(nM,1,'double');
REC_OSPB_dry      = nan(nM,1,'double');
G_ISPBwet_OSPBdry = nan(nM,1,'double');

for i = 1:nM
  fprintf('[Evaluate] Model %d/%d: %s\n', i, nM, S.names{i});

  x    = gdif_cache{i};                 % ROI-vector
  dm_i = get_model_margin(cfg, S, i);
  margin_used(i) = dm_i;

  prd = x(M) >= dm_i;                   % M-space predictions
  yr  = y_raw_vec(M);                   % M-space truth

  m_raw = compute_metrics(prd, yr);

  RAW_G(i)   = sqrt(max(m_raw.REC,0) * max(m_raw.SPEC,0));

  TP(i)=m_raw.TP; FP(i)=m_raw.FP; TN(i)=m_raw.TN; FN(i)=m_raw.FN;
  ACC(i)=m_raw.ACC; PR_r(i)=m_raw.PREC; RC(i)=m_raw.REC; F1(i)=m_raw.F1;
  TPFP(i)=m_raw.TP_FP; TNFN(i)=m_raw.TN_FN;

  % Regional recalls in ROI-vector space
  pred_full = false(numel(idx),1);
  pred_full(M) = prd;

  if N_wet_ISPB > 0
      TP_wet = sum(pred_full(mask_wet_ISPB_vec));
      REC_ISPB_wet(i) = TP_wet / max(N_wet_ISPB, eps);
  end

  if N_dry_OSPB > 0
      TN_dry = sum(~pred_full(mask_dry_OSPB_vec));
      REC_OSPB_dry(i) = TN_dry / max(N_dry_OSPB, eps);
  end

  if isfinite(REC_ISPB_wet(i)) && isfinite(REC_OSPB_dry(i))
      G_ISPBwet_OSPBdry(i) = sqrt(max(REC_ISPB_wet(i),0) * max(REC_OSPB_dry(i),0));
  end
end

if ~isfield(S,'titles') || numel(S.titles)~=nM, S.titles = S.names; end

results_table = table( ...
  S.titles(:), nnz(M)*ones(nM,1), ...
  TP, FP, TN, FN, ...
  ACC, PR_r, RC, F1, ...
  RAW_G, ...
  REC_ISPB_wet, REC_OSPB_dry, G_ISPBwet_OSPBdry, ...
  TPFP, TNFN, ...
  'VariableNames', {'Model','N_roi', ...
  'TP_raw','FP_raw','TN_raw','FN_raw', ...
  'ACC_raw','PREC_raw','REC_raw','F1_raw', ...
  'G_raw', ...
  'REC_ISPB_wet_raw','REC_OSPB_dry_raw','G_ISPBwet_OSPBdry_raw', ...
  'TP_FP_raw','TN_FN_raw'});

% Save point table
timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
out_csv = fullfile(cfg.outdir, sprintf('results_table_%s.csv', timestamp));
writetable(results_table, out_csv);
fprintf('Saved results table: %s\n', out_csv);

% -------------------- Expected metrics + bootstrap CIs --------------------
fprintf('Analytic expected metrics with %s bootstrap (probabilistic)...\n', cfg.uncertainty.bootstrap_mode);

sigG_M = STATIC.sigG_M;

unc_source = strings(nM,1);

stat_list = {'SPEC','REC','PREC','ACC','F1'};
EXP_raw = struct(); CI_raw = struct(); CIi_raw = struct();
for s = stat_list
    EXP_raw.(s{1}) = nan(nM,1);
    CI_raw.(s{1})  = nan(nM,2);
    CIi_raw.(s{1}) = nan(nM,2);
end
EXP_raw.G = nan(nM,1); CI_raw.G = nan(nM,2); CIi_raw.G = nan(nM,2);

EXP_raw.REC_ISPB_wet      = nan(nM,1);
EXP_raw.REC_OSPB_dry      = nan(nM,1);
EXP_raw.G_ISPBwet_OSPBdry = nan(nM,1);

CI_raw.REC_ISPB_wet      = nan(nM,2);
CI_raw.REC_OSPB_dry      = nan(nM,2);
CI_raw.G_ISPBwet_OSPBdry = nan(nM,2);

CIi_raw.REC_ISPB_wet      = nan(nM,2);
CIi_raw.REC_OSPB_dry      = nan(nM,2);
CIi_raw.G_ISPBwet_OSPBdry = nan(nM,2);

use_comp_boot = strcmpi(cfg.uncertainty.bootstrap_mode,'component') && ...
                isfield(S,'boot_groups') && S.boot_groups.nC>0;
G_groups = {};
if use_comp_boot, G_groups = S.boot_groups.groups; end

if isfield(cfg.uncertainty,'seed') && ~isempty(cfg.uncertainty.seed), rng(cfg.uncertainty.seed); end

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};

    gdif_M_all = double(gdif_cache{i}(M)); % M-space gdif

    % σ(model) in M-space
    sigM_M = zeros(NM,1,'like',gdif_M_all);

    if isfield(S,'unc') && isfield(S.unc, fld) && ~isempty(S.unc.(fld)) && isequal(size(S.unc.(fld)), size(S.H))
        tmp      = S.unc.(fld);
        sigM_roi = tmp(idx);
        sigM_M   = sigM_roi(M);
        unc_source(i) = "per-pixel";
    elseif strcmp(name_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        half_roi = 0.5*abs(S.bounds.Losing_max(idx) - S.bounds.Losing_min(idx));
        sigM_M   = (half_roi(M)) / 1.96;
        unc_source(i) = "bounds->sigma";
    else
        sigM_M(:) = STATIC.scalar_sigma_per_model(i);
        unc_source(i) = "scalar-fallback";
    end

    dm_i = get_model_margin(cfg, S, i);

    % Wet probability in M-space
    pw_all = wet_prob_with_margin(gdif_M_all, sigM_M, sigG_M, dm_i);

    V = isfinite(pw_all) & isfinite(y_raw_M);
    pw = pw_all(V);
    yr = y_raw_M(V);
    if isempty(pw), continue; end

    mask_wet_ISPB_valid = mask_wet_ISPB_M(V);
    mask_dry_OSPB_valid = mask_dry_OSPB_M(V);

    % Analytic expectations
    ETP  = sum(pw .* yr);          EFP  = sum(pw .* (1-yr));
    EFN  = sum((1-pw) .* yr);      ETN  = sum((1-pw) .* (1-yr));
    clamp = @(x) max(x, 0);
    [ETP,EFP,ETN,EFN] = deal(clamp(ETP), clamp(EFP), clamp(ETN), clamp(EFN));

    PREC_E = ETP / max(ETP + EFP, eps);
    SPEC_E = ETN / max(ETN + EFP, eps);
    REC_E  = ETP / max(ETP + EFN, eps);
    ACC_E  = (ETP + ETN) / max(numel(pw), eps);
    F1_E   = 2 * PREC_E * REC_E / max(PREC_E + REC_E, eps);

    EXP_raw.PREC(i)=PREC_E;
    EXP_raw.SPEC(i)=SPEC_E;
    EXP_raw.REC(i) =REC_E;
    EXP_raw.ACC(i) =ACC_E;
    EXP_raw.F1(i)  =F1_E;
    EXP_raw.G(i)   = sqrt(max(REC_E,0)*max(SPEC_E,0));

    % Regional expectations
    denom_wet_ISPB = sum( yr(mask_wet_ISPB_valid) );
    if denom_wet_ISPB > 0
        ETP_wet = sum( pw(mask_wet_ISPB_valid) .* yr(mask_wet_ISPB_valid) );
        EXP_raw.REC_ISPB_wet(i) = ETP_wet / max(denom_wet_ISPB, eps);
    end

    denom_dry_OSPB = sum( (1-yr(mask_dry_OSPB_valid)) );
    if denom_dry_OSPB > 0
        ETN_dry = sum( (1-pw(mask_dry_OSPB_valid)) .* (1-yr(mask_dry_OSPB_valid)) );
        EXP_raw.REC_OSPB_dry(i) = ETN_dry / max(denom_dry_OSPB, eps);
    end

    if isfinite(EXP_raw.REC_ISPB_wet(i)) && isfinite(EXP_raw.REC_OSPB_dry(i))
        EXP_raw.G_ISPBwet_OSPBdry(i) = sqrt(max(EXP_raw.REC_ISPB_wet(i),0) * max(EXP_raw.REC_OSPB_dry(i),0));
    end

    % Bootstrap replicates
    B  = cfg.uncertainty.n_boot;
    BR = struct();
    for s = stat_list, BR.(s{1}) = zeros(B,1); end
    BR.G = zeros(B,1);
    BR.REC_ISPB_wet      = nan(B,1);
    BR.REC_OSPB_dry      = nan(B,1);
    BR.G_ISPBwet_OSPBdry = nan(B,1);

    pred_mode = lower(string(cfg.uncertainty.pred_mode));
    use_noisy = strcmp(pred_mode,"noisy_threshold") && isfield(cfg.uncertainty,'noise_enable') && cfg.uncertainty.noise_enable;
    s_all = [];

    if strcmpi(cfg.uncertainty.bootstrap_mode,'pixel')

        for b = 1:B
            rb = randi(numel(yr), numel(yr), 1);
            y_b = yr(rb);

            if use_noisy
                if isempty(s_all), s_all = sqrt(max(sigM_M(V).^2 + sigG_M(V).^2, eps)); end
                epsi = randn(numel(rb),1) .* (cfg.uncertainty.noise_K .* s_all(rb));
                gjit = gdif_M_all(V); gjit = gjit(rb) + epsi;
                prd  = (gjit >= dm_i);
            else
                prd  = rand(numel(rb),1) < pw(rb);
            end

            mB = compute_metrics(prd, y_b);
            BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
            BR.G(b)   = sqrt(max(mB.REC,0)*max(mB.SPEC,0));

            mask_wet_b = mask_wet_ISPB_valid(rb);
            mask_dry_b = mask_dry_OSPB_valid(rb);

            denom_wet_b = sum( y_b & mask_wet_b );
            if denom_wet_b > 0
                TP_wet_b = sum( prd & y_b & mask_wet_b );
                BR.REC_ISPB_wet(b) = TP_wet_b / max(denom_wet_b, eps);
            end

            denom_dry_b = sum( ~y_b & mask_dry_b );
            if denom_dry_b > 0
                TN_dry_b = sum( ~prd & ~y_b & mask_dry_b );
                BR.REC_OSPB_dry(b) = TN_dry_b / max(denom_dry_b, eps);
            end

            if isfinite(BR.REC_ISPB_wet(b)) && isfinite(BR.REC_OSPB_dry(b))
                BR.G_ISPBwet_OSPBdry(b) = sqrt(max(BR.REC_ISPB_wet(b),0) * max(BR.REC_OSPB_dry(b),0));
            end
        end

    else
        % Component bootstrap: build V-indices per component (robust mapping)
        GV = {};
        if ~isempty(G_groups)

            NM_local = nnz(M);
            NV_local = nnz(V);

            roi_to_Midx = zeros(numel(M),1,'uint32');    % numel(M) == numel(idx)
            roi_to_Midx(M) = uint32(1:NM_local);

            Midx_to_Vidx = zeros(NM_local,1,'uint32');
            Midx_to_Vidx(V) = uint32(1:NV_local);

            GV = cell(size(G_groups));
            for cc = 1:numel(G_groups)
                roi_pos = G_groups{cc};                     % ROI indices
                midx    = double(roi_to_Midx(roi_pos));     % M indices (0 if not in M)
                midx    = midx(midx > 0);

                vidx    = double(Midx_to_Vidx(midx));       % V indices (0 if not in V)
                vidx    = vidx(vidx > 0);

                if ~isempty(vidx)
                    GV{cc} = vidx(:);
                end
            end
            GV = GV(~cellfun(@isempty,GV));
        end
        nC = numel(GV);

        if nC==0
            warning('[bootstrap] No valid components -> pixel bootstrap fallback (within model).');
            for b = 1:B
                rb  = randi(numel(yr), numel(yr), 1);
                y_b = yr(rb);

                if use_noisy
                    if isempty(s_all), s_all = sqrt(max(sigM_M(V).^2 + sigG_M(V).^2, eps)); end
                    epsi = randn(numel(rb),1) .* (cfg.uncertainty.noise_K .* s_all(rb));
                    gjit = gdif_M_all(V); gjit = gjit(rb) + epsi;
                    prd  = (gjit >= dm_i);
                else
                    prd  = rand(numel(rb),1) < pw(rb);
                end

                mB = compute_metrics(prd, y_b);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0) * max(mB.SPEC,0));
            end
        else
            mfrac = max(0, min(1, cfg.uncertainty.mfrac));

            if cfg.uncertainty.pps_components
                w = cellfun(@numel, GV); sW = sum(w);
                if sW>0, w = w / sW; else, w = []; end
            else
                w = [];
            end

            for b = 1:B
                mC = max(1, round(mfrac * nC));
                if isempty(w)
                    rbC = randi(nC, mC, 1);
                else
                    u = rand(mC,1); cdf = cumsum(w(:));
                    rbC = arrayfun(@(x) find(cdf>=x,1,'first'), u);
                end

                pickV = vertcat(GV{rbC});    % V indices directly
                pv = pw(pickV);
                yv = yr(pickV);

                if use_noisy
                    if isempty(s_all), s_all = sqrt(max(sigM_M(V).^2 + sigG_M(V).^2, eps)); end
                    epsi = randn(numel(pickV),1) .* (cfg.uncertainty.noise_K .* s_all(pickV));
                    gjit = gdif_M_all(V); gjit = gjit(pickV) + epsi;
                    prd  = (gjit >= dm_i);
                else
                    prd  = rand(numel(pv),1) < pv;
                end

                mB = compute_metrics(prd, yv);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0) * max(mB.SPEC,0));

                mask_wet_b = mask_wet_ISPB_valid(pickV);
                mask_dry_b = mask_dry_OSPB_valid(pickV);

                denom_wet_b = sum( yv & mask_wet_b );
                if denom_wet_b > 0
                    TP_wet_b = sum( prd & yv & mask_wet_b );
                    BR.REC_ISPB_wet(b) = TP_wet_b / max(denom_wet_b, eps);
                end

                denom_dry_b = sum( ~yv & mask_dry_b );
                if denom_dry_b > 0
                    TN_dry_b = sum( ~prd & ~yv & mask_dry_b );
                    BR.REC_OSPB_dry(b) = TN_dry_b / max(denom_dry_b, eps);
                end

                if isfinite(BR.REC_ISPB_wet(b)) && isfinite(BR.REC_OSPB_dry(b))
                    BR.G_ISPBwet_OSPBdry(b) = sqrt(max(BR.REC_ISPB_wet(b),0) * max(BR.REC_OSPB_dry(b),0));
                end
            end
        end
    end

    % CI summaries
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

% Build Gstats struct (per-run)
Gstats = struct();
for s = stat_list
    Gstats.([s{1} '_raw']) = struct('EXP', EXP_raw.(s{1}), 'CI', CI_raw.(s{1}), 'CIi', CIi_raw.(s{1}));
end
Gstats.G_raw = struct('EXP', EXP_raw.G, 'CI', CI_raw.G, 'CIi', CIi_raw.G);

Gstats.REC_ISPB_wet_raw = struct('EXP', EXP_raw.REC_ISPB_wet, 'CI', CI_raw.REC_ISPB_wet, 'CIi', CIi_raw.REC_ISPB_wet);
Gstats.REC_OSPB_dry_raw = struct('EXP', EXP_raw.REC_OSPB_dry, 'CI', CI_raw.REC_OSPB_dry, 'CIi', CIi_raw.REC_OSPB_dry);
Gstats.G_ISPBwet_OSPBdry_raw = struct('EXP', EXP_raw.G_ISPBwet_OSPBdry, 'CI', CI_raw.G_ISPBwet_OSPBdry, 'CIi', CIi_raw.G_ISPBwet_OSPBdry);

% Attach sigma summaries and save
results_table.sigmaGdel_mean   = S.sigma_gdif_mean;
results_table.sigmaGdel_median = S.sigma_gdif_median;
timestamp_sigma = char(datetime('now','Format','yyyyMMdd_HHmmss'));
out_csv_sigma = fullfile(cfg.outdir, sprintf('sigmaGdel_summary_%s.csv', timestamp_sigma));
writetable(results_table(:, {'Model','sigmaGdel_mean','sigmaGdel_median'}), out_csv_sigma);
fprintf('Saved sigma(GΔ) summary: %s\n', out_csv_sigma);

% CI summary CSV (compact)
try
    Name  = string(S.names(:));
    Title = string(S.titles(:));

    Tsum = table(Name, Title, ...
        i_pad(Gstats.SPEC_raw.EXP(:),nM), i_pad(Gstats.SPEC_raw.CI(:,1),nM), i_pad(Gstats.SPEC_raw.CI(:,2),nM), ...
        i_pad(Gstats.REC_raw.EXP(:),nM),  i_pad(Gstats.REC_raw.CI(:,1),nM),  i_pad(Gstats.REC_raw.CI(:,2),nM),  ...
        i_pad(Gstats.PREC_raw.EXP(:),nM), i_pad(Gstats.PREC_raw.CI(:,1),nM), i_pad(Gstats.PREC_raw.CI(:,2),nM), ...
        i_pad(Gstats.ACC_raw.EXP(:),nM),  i_pad(Gstats.ACC_raw.CI(:,1),nM),  i_pad(Gstats.ACC_raw.CI(:,2),nM),  ...
        i_pad(Gstats.F1_raw.EXP(:),nM),   i_pad(Gstats.F1_raw.CI(:,1),nM),   i_pad(Gstats.F1_raw.CI(:,2),nM),   ...
        'VariableNames', {'Name','Title', ...
          'SPEC_raw','SPEC_raw_lo','SPEC_raw_hi', ...
          'REC_raw','REC_raw_lo','REC_raw_hi', ...
          'PREC_raw','PREC_raw_lo','PREC_raw_hi', ...
          'ACC_raw','ACC_raw_lo','ACC_raw_hi', ...
          'F1_raw','F1_raw_lo','F1_raw_hi'});

    Tsum.G_raw    = i_pad(Gstats.G_raw.EXP(:),nM);
    Tsum.G_raw_lo = i_pad(Gstats.G_raw.CI(:,1),nM);
    Tsum.G_raw_hi = i_pad(Gstats.G_raw.CI(:,2),nM);

    Tsum.REC_ISPB_wet_raw    = i_pad(Gstats.REC_ISPB_wet_raw.EXP(:),nM);
    Tsum.REC_ISPB_wet_raw_lo = i_pad(Gstats.REC_ISPB_wet_raw.CI(:,1),nM);
    Tsum.REC_ISPB_wet_raw_hi = i_pad(Gstats.REC_ISPB_wet_raw.CI(:,2),nM);

    Tsum.REC_OSPB_dry_raw    = i_pad(Gstats.REC_OSPB_dry_raw.EXP(:),nM);
    Tsum.REC_OSPB_dry_raw_lo = i_pad(Gstats.REC_OSPB_dry_raw.CI(:,1),nM);
    Tsum.REC_OSPB_dry_raw_hi = i_pad(Gstats.REC_OSPB_dry_raw.CI(:,2),nM);

    Tsum.G_ISPBwet_OSPBdry_raw    = i_pad(Gstats.G_ISPBwet_OSPBdry_raw.EXP(:),nM);
    Tsum.G_ISPBwet_OSPBdry_raw_lo = i_pad(Gstats.G_ISPBwet_OSPBdry_raw.CI(:,1),nM);
    Tsum.G_ISPBwet_OSPBdry_raw_hi = i_pad(Gstats.G_ISPBwet_OSPBdry_raw.CI(:,2),nM);

    timestamp_ci = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    csv_path = fullfile(cfg.outdir, sprintf('ci_summary_%s.csv', timestamp_ci));
    writetable(Tsum, csv_path);
    fprintf('Saved CI summary CSV: %s\n', csv_path);
catch ME
    fprintf('[warn] CI summary CSV failed (non-fatal): %s\n', ME.message);
end

% -------------------- Optional lightweight artifacts (per run) --------------------
if cfg.save_artifacts
    save_lightweight_plot_bundle(cfg, S, y_raw_full, Gstats);
end

% -------------------- Return struct (per run) --------------------
OUT = struct();
OUT.cfg = cfg;
OUT.spec_thresh = cfg.spec_thresh;
OUT.results_table = results_table;
OUT.Gstats = Gstats;
OUT.margin_used = margin_used;
OUT.unc_source = unc_source;

end

%% ========================================================================
%  ARTIFACT BUNDLE (KEEP SMALL; PER RUN)
% ========================================================================
function save_lightweight_plot_bundle(cfg, S, y_raw_full, Gstats)

cfg.artifact_dir = fullfile(cfg.outdir,'artifacts');
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end

ds = 8;
r = 1:ds:size(S.H,1);
c = 1:ds:size(S.H,2);

S_plot = struct();
S_plot.includes = {'X_km','Y_km','H','viz_mask','sink_mask_comp','comp_id', ...
                   'roi_mask','y_raw_full','mask_ISPB','mask_OSPB', ...
                   'mask_SPB','ds_stride','names','titles','Gstats', ...
                   'model_cvals_median','model_cvals','model_quantiles'};

S_plot.names  = S.names;
if isfield(S,'titles') && numel(S.titles)==numel(S.names), S_plot.titles = S.titles; else, S_plot.titles = S.names; end
S_plot.run_id = cfg.run_id;
S_plot.Gstats = Gstats;

if isfield(S,'model_cvals_median'), S_plot.model_cvals_median = S.model_cvals_median(:); end
if isfield(S,'model_cvals'),        S_plot.model_cvals        = S.model_cvals(:);        end
if isfield(S,'model_quantiles'),    S_plot.model_quantiles    = S.model_quantiles;       end

S_plot.X_km     = single(S.X_km(r,c));
S_plot.Y_km     = single(S.Y_km(r,c));
S_plot.H        = single(S.H(r,c));
S_plot.viz_mask = logical(S.viz_mask(r,c));
S_plot.sink_mask_comp = logical(S.sink_mask_comp(r,c));
S_plot.comp_id        = uint32(S.comp_id(r,c));
S_plot.roi_mask       = logical(S.roi_mask(r,c));

if ~isempty(y_raw_full)
    S_plot.y_raw_full = logical(y_raw_full(r,c));
end

if isfield(S,'mask_ISPB'), S_plot.mask_ISPB = logical(S.mask_ISPB(r,c)); end
if isfield(S,'mask_OSPB'), S_plot.mask_OSPB = logical(S.mask_OSPB(r,c)); end
if isfield(S,'mask_SPB'),  S_plot.mask_SPB  = logical(S.mask_SPB(r,c));  end

S_plot.ds_stride = uint16(ds);

cfg_plot = struct();
cfg_plot.outdir = cfg.outdir;
cfg_plot.run_id = cfg.run_id;
cfg_plot.font_size = cfg.font_size;
cfg_plot.overlay_alpha = cfg.overlay_alpha;

stamp  = char(datetime('now','Format','yyyyMMdd_HHmmss'));
f_plot = fullfile(cfg.artifact_dir, sprintf('S_plot_small_%s.mat', stamp));
save(f_plot, 'S_plot', 'cfg_plot', '-v7');
fprintf('[artifacts] Saved lightweight plot file: %s\n', f_plot);

end

%% ========================================================================
%  MASK STITCH / RESIZE HELPERS
% ========================================================================
function S = stitch_spb_masks_to_grid(S, m)
% Reconstruct cropped SPB masks onto full S.Xgrid size (ONLY for regional metrics)
sz = size(S.Xgrid);

S.mask_ISPB = false(sz);
S.mask_OSPB = false(sz);
S.mask_SPB  = false(sz);

required = {'rmin','rmax','cmin','cmax','mask_ISPB_c','mask_OSPB_c','mask_SPB_c'};
miss = required(~isfield(m,required));
if ~isempty(miss)
    error('spb_masks.mat missing fields: %s', strjoin(miss,', '));
end

S.mask_ISPB(m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_ISPB_c));
S.mask_OSPB(m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_OSPB_c));
S.mask_SPB( m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_SPB_c));

% REG_MASK removed -> rect_mask not needed; keep for compatibility but all-true.
S.rect_mask = true(sz);

end

function [mask_out, did] = force_mask_to_size(mask_in, target_sz, name)
did = false;
mask_out = mask_in;
if isempty(mask_in)
    mask_out = false(target_sz);
    did = true;
    warning('%s was empty; created false mask.', name);
    return
end
if isequal(size(mask_in), target_sz)
    return
end
if isequal(size(mask_in), fliplr(target_sz))
    mask_in = mask_in.';
end
if ~isequal(size(mask_in), target_sz)
    did = true;
    warning('%s size %s != target %s -> resizing (nearest).', name, mat2str(size(mask_in)), mat2str(target_sz));
    mask_out = logical(imresize(logical(mask_in), target_sz, 'nearest'));
else
    mask_out = logical(mask_in);
end
end

%% ========================================================================
%  DEPENDENCY / PATH HELPERS
% ========================================================================
function check_dependencies_local(opts)
% Best-effort check. If you have your own check_dependencies.m, it will be used.
if exist('check_dependencies','file') == 2
    check_dependencies(opts);
    return
end

if ~isfield(opts,'verbose'), opts.verbose = true; end

hasTB = @(names) any(cellfun(@(n) license('test', n), names));
hasImage = hasTB({'Image_Toolbox','image_toolbox'});
hasStats = hasTB({'Statistics_Toolbox','statistics_toolbox'});
hasPCT   = hasTB({'Distrib_Computing_Toolbox','distrib_computing_toolbox'});

if opts.verbose
    if hasImage
        fprintf('[deps] OK: Image Processing Toolbox\n');
    else
        warning('[deps] Missing Image Processing Toolbox (needed for bwlabel/imresize).');
    end

    if hasStats
        fprintf('[deps] OK: Statistics and Machine Learning Toolbox\n');
    else
        warning('[deps] Missing Statistics and Machine Learning Toolbox (needed for prctile).');
    end

    if hasPCT
        fprintf('[deps] OK (optional): Parallel Computing Toolbox\n');
    else
        fprintf('[deps] Parallel Computing Toolbox not found (optional).\n');
    end
end
end

function require_on_path(fn_list)
missing = {};
for i = 1:numel(fn_list)
    if exist(fn_list{i},'file') ~= 2
        missing{end+1} = fn_list{i}; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error(['Missing required function(s): ' strjoin(missing,', ') ...
           '. Add your project code folders to the MATLAB path.']);
end
end

%% ========================================================================
%  SMALL HELPERS
% ========================================================================
function m = compute_metrics(pred, truth)
pred  = logical(pred(:)); truth = logical(truth(:));
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
m = struct('N',N, 'TP',TP,'FP',FP,'TN',TN,'FN',FN, ...
           'ACC',ACC,'PREC',PREC,'SPEC',SPEC,'REC',REC,'F1',F1, ...
           'TP_FP', TP / max(FP, eps), ...
           'TN_FN', TN / max(FN, eps));
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

haveNames = isfield(S,'names') && iscell(S.names) && ~isempty(S.names);
if haveNames
    is_flat = startsWith(S.names,'Flat_'); is_flat = is_flat(:);
    nf_mask = ~is_flat; if ~any(nf_mask), nf_mask = true(size(is_flat)); end
else
    nf_mask = true(1,1);
end

haveMean   = isfield(S,'sigma_gdif_mean')   && ~isempty(S.sigma_gdif_mean);
haveMedian = isfield(S,'sigma_gdif_median') && ~isempty(S.sigma_gdif_median);
safe_med = @(v,mask) median(v(mask), 'omitnan');

switch mode
    case "auto_mean"
        if isfield(S,'suggested_margin_mean') && isfinite(S.suggested_margin_mean)
            dm = S.suggested_margin_mean;
        elseif haveMean
            dm = safe_med(S.sigma_gdif_mean(:), nf_mask);
        else
            dm = cfg.decision_margin;
        end
    case "auto_median"
        if isfield(S,'suggested_margin_median') && isfinite(S.suggested_margin_median)
            dm = S.suggested_margin_median;
        elseif haveMedian
            dm = safe_med(S.sigma_gdif_median(:), nf_mask);
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

function v = i_pad(v, n)
v = v(:);
if numel(v) < n, v(end+1:n,1) = NaN; end
if numel(v) > n, v = v(1:n); end
end
