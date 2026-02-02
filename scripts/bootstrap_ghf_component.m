%% ================================================================
%  GHF model evaluation 
%  ================================================================
function S = bootstrap_ghf_component(cfg_override)
check_dependencies(struct('autoInstall',true,'needAMT',false,'verbose',true));
addpath('datasets');
dbstop if error
run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));  

%% -------------------- CONFIG --------------------
cfg = struct( ...
  'spec_thresh',        0.2, ...   % Specularity content threshold
  'v_keep',             10,  ...   % Ice velocity threshold (m/yr)
  'decision_margin',    0,   ...   % Fixed margin: wet if (GΔ >= 0)
  'decision_margin_mode','fixed', ... % 'fixed' | 'auto_mean' | 'auto_median' | 'auto_per_model_mean' | 'auto_per_model_median'
  'outdir',             'figs_out', ...
  'overwrite',          true, ...
  'overlay_alpha',      0.4, ...
  'font_size',          20, ...
  'to_single',          true, ...
  'sweep_mode',         true ...
);

cfg.save_artifacts = true;
cfg.artifact_dir   = fullfile(cfg.outdir,'artifacts');
if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end

cfg.run_id = run_id;

% === Region split from boundary polyline (ISPB/OSPB/all) ===
cfg.region = struct( ...
  'mode',          'ALL', ...         % 'ISPB' | 'OSPB' | 'SPB' | 'ALL'
  'left_label',    'ISPB', ...
  'right_label',   'OSPB', ...
  'cache',         true ...
);

% Sinks
cfg.sinks = struct('marker','.', 'size',4, 'color',[0.05 0.05 0.05], 'alpha',0.85);

cfg.skip_static_build = true;

% ---------------- Uncertainty — analytic expectation + sink bootstrap -----
cfg.uncertainty = struct( ...
  'mode','analytic', ...
  'bootstrap_mode','component', ...    % 'pixel' | 'component'
  'n_boot',10, ...
  'seed',42, ...
  'mfrac',0.7, ...
  'band_main',0.95, ...
  'band_inner',0.50, ...
  'noise_enable',false, ...
  'noise_K',1, ...
  'pred_mode','bernoulli', ...     % 'bernoulli' | 'noisy_threshold'
  'class_balance',false, ...       % pixel bootstrap: stratified by class
  'pps_components',false ...       % component bootstrap: PPS by component size
);

cfg.plot_labels = struct('avoid_overlap',true,'max_iter',200,'step',0.01,'pad',0.005);
cfg.synthetic_only = false;

if nargin > 0 && ~isempty(cfg_override)
    cfg = merge_structs(cfg, cfg_override);
end
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end

% =========================== Load data ===========================
fprintf('Loading data ...\n');
datafile = '/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/gmin_data.mat';

if exist(datafile,'file') == 2
    fprintf('[bootstrap] Loading %s ... ', datafile);
    L = load(datafile,'S');
    S = L.S;
    clear L;
    fprintf('ok.\n');
else
    warning('Data file not found. Building now...');

    S = build_data_core('to_single',true,'save',true, ...
                    'outdir','datasets','outfile','gmin_data.mat');
    S = add_g_fields(fullfile('datasets','gmin_data.mat'), ...
                 'advection_file','datasets/longitudinal_advection_maps.mat', ...
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
fprintf('[bootstrap] Grid: %dx%d | models: %d\n', size(S.Xgrid,1), size(S.Xgrid,2), numel(S.names));

% Finalize model list
nM = numel(S.names); S.model_cvals = nan(nM,1);

%% =========================== Load SPB masks ===========================
script_dir  = fileparts(mfilename('fullpath'));
dataset_dir = fullfile(script_dir, 'datasets');
if ~exist(dataset_dir,'dir'), mkdir(dataset_dir); end

datafile = fullfile(dataset_dir, 'spb_masks.mat');

if exist(datafile,'file') == 2
    m = load(datafile);
else
    old = pwd;
    cd(script_dir);          % so build_spb_masks() saves to script_dir/datasets/...
    build_spb_masks();
    cd(old);

    if exist(datafile,'file') ~= 2
        error('build_spb_masks() did not create expected file: %s', datafile);
    end
    m = load(datafile);
end

S.mask_ISPB = false(size(S.Xgrid));
S.mask_OSPB = false(size(S.Xgrid));
S.mask_SPB  = false(size(S.Xgrid));

S.mask_ISPB(m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_ISPB_c));
S.mask_OSPB(m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_OSPB_c));
S.mask_SPB( m.rmin:m.rmax, m.cmin:m.cmax) = logical(full(m.mask_SPB_c));
S.rect_mask = logical(m.rect_mask);

% switch lower(cfg.region.mode)
%   case 'ispb',  REG_MASK = S.mask_ISPB & S.rect_mask;
%   case 'ospb',  REG_MASK = S.mask_OSPB & S.rect_mask;
%   case 'spb',   REG_MASK = S.mask_SPB  & S.rect_mask;
%   case 'all',   REG_MASK = S.rect_mask;
%   otherwise,    error('cfg.region.mode must be ISPB | OSPB | SPB | ALL');
% end

%% ========== ROI masks (sinks + slow flow + region) ==========
if exist('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_new.mat','file')==2 && exist('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat','file')==2
    t = load('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_new.mat', 'sink_mask'); S.sink_mask = t.sink_mask; clear t;    
    t = load('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat', 'comp_id','sink_mask_comp'); 
    S.comp_id = t.comp_id; S.sink_mask_comp = t.sink_mask_comp; clear t;
else
    make_sink_mask();
    t = load('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_new.mat', 'sink_mask'); S.sink_mask = t.sink_mask; clear t;  
    t = load('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat', 'comp_id','sink_mask_comp'); 
    S.comp_id = t.comp_id; S.sink_mask_comp = t.sink_mask_comp; clear t;
end

valid_mask = isfinite(S.Q) & ~S.spec_invalid & isfinite(S.H) & isfinite(S.icevel);
slow_mask  = valid_mask & ~(S.icevel > cfg.v_keep);

%% --- Build REG_MASK on the canonical grid (self-contained) ---
% Default: full domain
REG_MASK = true(size(S.H));

% OPTION A (recommended): bounding box in km: [xmin xmax ymin ymax]
% Example: cfg.reg_bbox_km = [-50 250 -120 160];  % <-- YOU set this
if isfield(cfg,'reg_bbox_km') && ~isempty(cfg.reg_bbox_km)
    b = cfg.reg_bbox_km(:).';
    assert(numel(b)==4, 'cfg.reg_bbox_km must be [xmin xmax ymin ymax] in km');
    xmin = b(1)*1000; xmax = b(2)*1000;
    ymin = b(3)*1000; ymax = b(4)*1000;

    REG_MASK = (S.Xgrid >= xmin & S.Xgrid <= xmax & ...
                S.Ygrid >= ymin & S.Ygrid <= ymax);
end

% OPTION B: polygon vertices in km (overrides bbox if provided)
% Example:
% cfg.reg_poly_km = [x1 y1; x2 y2; ... ; xN yN];  % km
if isfield(cfg,'reg_poly_km') && ~isempty(cfg.reg_poly_km)
    P = cfg.reg_poly_km;
    assert(size(P,2)==2, 'cfg.reg_poly_km must be Nx2 [x_km y_km]');
    xv = P(:,1)*1000;  yv = P(:,2)*1000;
    REG_MASK = inpolygon(S.Xgrid, S.Ygrid, xv, yv);
end

% Safety: always logical and correct size
REG_MASK = logical(REG_MASK);
assert(isequal(size(REG_MASK), size(S.H)), 'REG_MASK size mismatch');

% --- Now safe ---
S.roi_mask = slow_mask & REG_MASK & S.sink_mask_comp;
 
S.viz_mask = slow_mask & REG_MASK;                      

% ROI-vector indexing
idx = find(S.roi_mask);
S.eval_mask_full = S.roi_mask;          
S.eval_mask      = true(numel(idx),1);  

% -- MEDIAN domain 
if isfield(S,'viz_mask') && ~isempty(S.viz_mask)
    med_mask = S.viz_mask;                      
else
    % Fallback if viz_mask wasn't built for some reason:
    valid_mask = isfinite(S.Q) & (S.Q ~= 9999) & isfinite(S.H) & isfinite(S.icevel);
    slow_mask  = valid_mask & ~(S.icevel > cfg.v_keep);
    med_mask   = slow_mask & REG_MASK;
end
idx_med = find(med_mask);

% 'Truth' labels
spec_roi  = S.Q(idx);
spec_ok   = isfinite(spec_roi) & (spec_roi ~= 9999) & (spec_roi > cfg.spec_thresh);

y_raw_vec        = spec_ok;
S.y_raw_vec      = y_raw_vec(:);
S.y_raw_full     = false(size(S.Q));
S.y_raw_full(idx)= S.y_raw_vec;

% Grouping for component bootstrap
S.boot_groups = struct('nC',0,'groups',{{}},'uniqC',[]);
if strcmpi(cfg.uncertainty.bootstrap_mode,'component')
    Mglob = S.eval_mask;
    compM = S.comp_id(idx); compM = compM(Mglob);
    if any(compM)
        [uniqC,~,gidx] = unique(compM(compM>0));
        groups = accumarray(gidx, find(compM>0), [], @(v){v});
        S.boot_groups.uniqC  = uniqC;
        S.boot_groups.groups = groups;
        S.boot_groups.nC     = numel(groups);
    end
end

validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);

if ~isfield(S,'comp_id') || isempty(S.comp_id)
    if exist('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat','file')
        tt = load('/disk/kea/WAIS/home/wais/users/mek/gmin_code/for_pub/datasets/sink_mask_comp.mat');
        if isfield(tt,'comp_id') && isequal(size(tt.comp_id), size(S.sink_mask_comp))
            S.comp_id = tt.comp_id;
        else
            S.comp_id = bwlabel(S.sink_mask_comp, 8);
        end
    else
        S.comp_id = bwlabel(S.sink_mask_comp, 8);
    end
end

%% >>> SANITY & COVERAGE
fprintf('\n---------- SANITY & COVERAGE ----------\n');
[nr,nc] = size(S.Xgrid);
fprintf('[grid] size=%dx%d  (%.2f Mpx)\n', nr, nc, numel(S.Xgrid)/1e6);

N_roi   = nnz(S.roi_mask);
N_viz   = nnz(S.viz_mask);
N_reg   = nnz(REG_MASK);
fprintf('[masks] REG=%d, VIZ=%d, ROI=%d  (ROI/REG=%.1f%%)\n', N_reg, N_viz, N_roi, 100*N_roi/max(N_reg,1));

roi_spec = S.Q(S.roi_mask);
roi_vel  = S.icevel(S.roi_mask);
roi_H    = S.H(S.roi_mask);

pct_nan_Q  = 100*mean(~isfinite(roi_spec));
pct_nan_v  = 100*mean(~isfinite(roi_vel));
pct_nan_H  = 100*mean(~isfinite(roi_H));
fprintf('[roi] NaN%% — Q: %.2f%%, vel: %.2f%%, H: %.2f%%\n', pct_nan_Q, pct_nan_v, pct_nan_H);

p_base = 100*mean(S.y_raw_vec);
fprintf('[labels] wet base rate in ROI-vector M-space: %.3f%% (N=%d)\n', p_base, numel(S.y_raw_vec));
fprintf('[flow] slow threshold v_keep=%.1f m/yr  | ROI is fully slow-flow by construction.\n', cfg.v_keep);

if isfield(S,'comp_id') && ~isempty(S.comp_id)
    comp_roi = S.comp_id(S.roi_mask);
    comp_roi = comp_roi(comp_roi>0);
    if ~isempty(comp_roi)
        [u,~,g] = unique(comp_roi);
        sizes = accumarray(g,1);
        q = quantile(sizes,[0 .25 .5 .75 .9 .95 1]);
        fprintf('[components] nC=%d | size{min, Q1, med, Q3, P90, P95, max} = [%d, %d, %d, %d, %d, %d, %d]\n', ...
            numel(u), q(1), q(2), q(3), q(4), q(5), q(6), q(7));
        if numel(u) < 50
            warning('[components] Low component count (%d). CI may be unstable; consider pixel bootstrap or looser gates.', numel(u));
        end
    else
        warning('[components] No component IDs > 0 inside ROI.');
    end
else
    warning('[components] comp_id missing; component bootstrap unavailable.');
end
fprintf('--------------------------------------\n\n');

%% ----------------- Per-model center values & GΔ cache -------------------
% reuse idx = find(S.roi_mask);
idx = find(S.roi_mask);

% Update model center values
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
        S.model_cvals_median(i) = v;          % flats get their scalar
        S.model_cvals(i)        = v;          % keep legacy field in sync
    else
        S.model_cvals_median(i) = median(Mi(idx_med), 'omitnan');  
        S.model_cvals(i)        = S.model_cvals_median(i);
    end
end
S.model_quantiles = struct('q50', S.model_cvals_median(:));

% GΔ cache in ROI-vector space
gdif_cache = cell(nM,1);
for i = 1:nM
    fld = validNames{i};
    gi  = S.Gdiff.(fld);
    gdif_cache{i} = gi(idx);
end

% ---- regional-class masks ----
mask_ISPB_vec = false(numel(idx),1);
mask_OSPB_vec = false(numel(idx),1);

if isfield(S,'mask_ISPB') && ~isempty(S.mask_ISPB)
    mask_ISPB_vec = S.mask_ISPB(idx);
end
if isfield(S,'mask_OSPB') && ~isempty(S.mask_OSPB)
    mask_OSPB_vec = S.mask_OSPB(idx);
end

wet_vec = logical(S.y_raw_vec(:));   % truth: wet = 1, dry = 0
dry_vec = ~wet_vec;

% Restrict to eval_mask
mask_wet_ISPB = S.eval_mask & mask_ISPB_vec & wet_vec;
mask_dry_OSPB = S.eval_mask & mask_OSPB_vec & dry_vec;

N_wet_ISPB = nnz(mask_wet_ISPB);
N_dry_OSPB = nnz(mask_dry_OSPB);

% Stash
S.mask_wet_ISPB_vec = mask_wet_ISPB;
S.mask_dry_OSPB_vec = mask_dry_OSPB;
S.N_wet_ISPB = N_wet_ISPB;
S.N_dry_OSPB = N_dry_OSPB;

%% -------------------- EVALUATE (point metrics, no uncertainty) ------------
fprintf('Evaluating point metrics on ROI...\n');

TP  = zeros(nM,1,'double'); FP  = TP; TN = TP; FN = TP;
ACC = TP; PR_r = TP; RC = TP; F1 = TP; TPFP=TP; TNFN=TP;
S.margin_used = nan(nM,1);

RAW_SPEC = nan(nM,1); RAW_REC = nan(nM,1); RAW_PREC = nan(nM,1);
RAW_ACC  = nan(nM,1); RAW_F1  = nan(nM,1); RAW_G = sqrt(max(RAW_SPEC,0) .* max(RAW_REC,0));

% regional-recall metric storage
REC_ISPB_wet      = nan(nM,1,'double');   % recall of wet in ISPB
REC_OSPB_dry      = nan(nM,1,'double');   % recall of dry in OSPB
G_ISPBwet_OSPBdry = nan(nM,1,'double');   % geometric mean of the two

for i = 1:nM
  fprintf('[Evaluate] Model %d/%d: %s\n', i, nM, S.names{i});
  x  = gdif_cache{i};
  Mv = S.eval_mask;
  dm_i = get_model_margin(cfg, S, i);
  S.margin_used(i) = dm_i;
  prd  = x(Mv) >= dm_i;      
  yr   = S.y_raw_vec(Mv);    

  m_raw = compute_metrics(prd, yr);

  RAW_SPEC(i)=m_raw.SPEC; RAW_REC(i)=m_raw.REC; RAW_PREC(i)=m_raw.PREC;
  RAW_ACC(i) =m_raw.ACC;  RAW_F1(i) =m_raw.F1;

  TP(i)=m_raw.TP; FP(i)=m_raw.FP; TN(i)=m_raw.TN; FN(i)=m_raw.FN;
  ACC(i)=m_raw.ACC; PR_r(i)=m_raw.PREC; RC(i)=m_raw.REC; F1(i)=m_raw.F1; 
  TPFP(i)=m_raw.TP_FP; TNFN(i)=m_raw.TN_FN;

  % ---- ISPB-wet / OSPB-dry metric ----
  pred_full = false(numel(idx),1);
  pred_full(Mv) = prd;

  if N_wet_ISPB > 0
      TP_wet_ISPB = sum(pred_full(mask_wet_ISPB));      
      REC_ISPB_wet(i) = TP_wet_ISPB / max(N_wet_ISPB, eps);
  else
      REC_ISPB_wet(i) = NaN;
  end

  if N_dry_OSPB > 0
      TN_dry_OSPB = sum(~pred_full(mask_dry_OSPB));     
      REC_OSPB_dry(i) = TN_dry_OSPB / max(N_dry_OSPB, eps);
  else
      REC_OSPB_dry(i) = NaN;
  end

  if isfinite(REC_ISPB_wet(i)) && isfinite(REC_OSPB_dry(i))
      G_ISPBwet_OSPBdry(i) = sqrt( max(REC_ISPB_wet(i),0) * max(REC_OSPB_dry(i),0) );
  else
      G_ISPBwet_OSPBdry(i) = NaN;
  end
end

if ~isfield(S,'titles') || numel(S.titles)~=nM, S.titles = S.names; end

S.results_table = table( ...
  S.titles(:), nnz(S.eval_mask)*ones(nM,1), ...
  TP, FP, TN, FN, ...
  ACC, PR_r, RC, F1, ...
  RAW_G, ...                               % overall G = sqrt(REC*SPEC)
  REC_ISPB_wet, REC_OSPB_dry, G_ISPBwet_OSPBdry, ...  
  TPFP, TNFN, ...
  'VariableNames', {'Model','N_roi', ...
  'TP_raw','FP_raw','TN_raw','FN_raw', ...
  'ACC_raw','PREC_raw','REC_raw','F1_raw', ...
  'G_raw', ...                                      % overall G
  'REC_ISPB_wet_raw','REC_OSPB_dry_raw','G_ISPBwet_OSPBdry_raw', ... 
  'TP_FP_raw','TN_FN_raw'});

% stash
S.REC_ISPB_wet      = REC_ISPB_wet;
S.REC_OSPB_dry      = REC_OSPB_dry;
S.G_ISPBwet_OSPBdry = G_ISPBwet_OSPBdry;

timestamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
out_csv = fullfile(cfg.outdir, sprintf('results_table%s_%s.csv', region_suffix(cfg), timestamp));
writetable(S.results_table, out_csv);
fprintf('Saved results table: %s\n', out_csv);

%% ---- EXPECTED METRICS + (PIXEL/COMPONENT) BOOTSTRAP ----------------------
fprintf('Analytic expected metrics with %s bootstrap (probabilistic)...\n', cfg.uncertainty.bootstrap_mode);

stat_list = {'SPEC','REC','PREC','ACC','F1'};

% Build once: ROI indexing and Gmin σ within ROI-vector M-space
M        = S.eval_mask;                 
y_raw_M  = double(S.y_raw_full(idx));   
y_raw_M  = y_raw_M(M);
NM       = nnz(M);

% regional masks in M-space 
mask_wet_ISPB_M = S.mask_wet_ISPB_vec(M);
mask_dry_OSPB_M = S.mask_dry_OSPB_vec(M);

% σ(Gmin) per pixel in ROI-vector 
if isfield(S,'sigma_gmin') && ~isempty(S.sigma_gmin) && isequal(size(S.sigma_gmin), size(S.H))
    sigG_M_full = S.sigma_gmin(idx);
    sigG_M      = sigG_M_full(M);
else
    sigG_M = zeros(NM,1,'like',y_raw_M);
end

% Prepare per-model caches and recompute σ(GΔ) aggregates inside ROI
validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);
S.sigma_gdif_mean   = nan(nM,1);
S.sigma_gdif_median = nan(nM,1);

scalar_sigma_per_model = nan(nM,1);   % final scalar fallback per model

if isfield(cfg.uncertainty,'seed') && ~isempty(cfg.uncertainty.seed), rng(cfg.uncertainty.seed); end

for i = 1:nM
    name_i = S.names{i};
    fld    = validNames{i};

    % Attempt per-pixel σ(model) first
    sigM_grid = [];
    if isfield(S,'unc') && isfield(S.unc,fld) && isequal(size(S.unc.(fld)), size(S.H))
        sigM_grid = S.unc.(fld);  % σ(model) on grid
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
        % No per-pixel map —> scalar σ(model)
        guess = NaN;
        if isfield(S,'sigma_model_scalar') && numel(S.sigma_model_scalar)>=i
            guess = S.sigma_model_scalar(i);
        elseif isfield(S,'sigma_gdif_mean') && numel(S.sigma_gdif_mean)>=i && isfinite(S.sigma_gdif_mean(i))
            guess = S.sigma_gdif_mean(i);
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

% existing global G
EXP_raw.G = nan(nM,1);
CI_raw.G  = nan(nM,2);
CIi_raw.G = nan(nM,2);

% regional metrics
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
        sigM_M = tmp(idx); 
        sigM_M = sigM_M(M);
    elseif strcmp(name_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        half_roi = 0.5*abs(S.bounds.Losing_max(idx) - S.bounds.Losing_min(idx));
        sigM_M   = (half_roi(M)) / 1.96;
    else
        sigM_M(:) = scalar_sigma_per_model(i);
    end

    if isfield(S,'unc') && isfield(S.unc, fld) && ~isempty(S.unc.(fld)) && isequal(size(S.unc.(fld)), size(S.H))
        S.unc_source(i) = "per-pixel";
    elseif strcmp(name_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        S.unc_source(i) = "bounds→sigma";
    else
        S.unc_source(i) = "scalar-fallback";
    end

    fprintf('[unc] sources: per-pixel=%d, bounds→sigma=%d, scalar-fallback=%d\n', ...
    nnz(S.unc_source=="per-pixel"), nnz(S.unc_source=="bounds→sigma"), nnz(S.unc_source=="scalar-fallback"));

    dm_i = get_model_margin(cfg, S, i);  
    pw_all = wet_prob_with_margin(gdif_M_all, sigM_M, sigG_M, dm_i);
    V      = isfinite(pw_all) & isfinite(y_raw_M);

    pw = pw_all(V); 
    yr = y_raw_M(V); 
    Nf = numel(pw); 
    if Nf==0, continue; end

    % valid masks
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
    ACC_E  = (ETP + ETN) / max(Nf, eps);
    F1_E   = 2 * PREC_E * REC_E / max(PREC_E + REC_E, eps);

    EXP_raw.PREC(i)=PREC_E; 
    EXP_raw.SPEC(i)=SPEC_E; 
    EXP_raw.REC(i) =REC_E; 
    EXP_raw.ACC(i) =ACC_E; 
    EXP_raw.F1(i)  =F1_E;

    % ---- regional expectations ----
    % expected ISPB wet recall 
    denom_wet_ISPB = sum( yr(mask_wet_ISPB_valid) );  
    if denom_wet_ISPB > 0
        ETP_wet_ISPB = sum( pw(mask_wet_ISPB_valid) .* yr(mask_wet_ISPB_valid) );
        REC_ISPB_wet_E = ETP_wet_ISPB / max(denom_wet_ISPB, eps);
    else
        REC_ISPB_wet_E = NaN;
    end

    % expected OSPB dry recall
    denom_dry_OSPB = sum( (1-yr(mask_dry_OSPB_valid)) );
    if denom_dry_OSPB > 0
        ETN_dry_OSPB = sum( (1-pw(mask_dry_OSPB_valid)) .* (1-yr(mask_dry_OSPB_valid)) );
        REC_OSPB_dry_E = ETN_dry_OSPB / max(denom_dry_OSPB, eps);
    else
        REC_OSPB_dry_E = NaN;
    end

    EXP_raw.REC_ISPB_wet(i) = REC_ISPB_wet_E;
    EXP_raw.REC_OSPB_dry(i) = REC_OSPB_dry_E;

    if isfinite(REC_ISPB_wet_E) && isfinite(REC_OSPB_dry_E)
        EXP_raw.G_ISPBwet_OSPBdry(i) = sqrt(max(REC_ISPB_wet_E,0) * max(REC_OSPB_dry_E,0));
    else
        EXP_raw.G_ISPBwet_OSPBdry(i) = NaN;
    end

    % ------------------- Bootstrap replicates -------------------
    B  = cfg.uncertainty.n_boot; 
    BR = struct(); 
    for s = stat_list, BR.(s{1}) = zeros(B,1); end
    BR.G = zeros(B,1); 
    BR.REC_ISPB_wet      = zeros(B,1);
    BR.REC_OSPB_dry      = zeros(B,1);
    BR.G_ISPBwet_OSPBdry = zeros(B,1);

    % Shared helpers for both modes
    pred_mode = lower(string(cfg.uncertainty.pred_mode));
    use_noisy = strcmp(pred_mode,"noisy_threshold") && isfield(cfg.uncertainty,'noise_enable') && cfg.uncertainty.noise_enable;
    s_all = [];  % lazily compute when needed

    if strcmpi(cfg.uncertainty.bootstrap_mode,'pixel')
        % -------- Pixel bootstrap ----------
        for b = 1:B
            pos = find(yr==1); 
            neg = find(yr==0);
            nPos = numel(pos); 
            nNeg = numel(neg);
            Nf   = nPos + nNeg;

            if cfg.uncertainty.class_balance && nPos>0 && nNeg>0
                rPos = pos(randi(nPos, nPos, 1));   
                rNeg = neg(randi(nNeg, nNeg, 1));   
                rb   = [rPos; rNeg];
                rb   = rb(randperm(Nf));
            else
                rb   = randi(Nf, Nf, 1);
            end

            y_b = yr(rb);

            if use_noisy
                if isempty(s_all), s_all = sqrt(max(sigM_M.^2 + sigG_M.^2, eps)); end
                epsi = randn(numel(rb),1) .* (cfg.uncertainty.noise_K .* s_all(rb));
                gjit = gdif_M_all(V); gjit = gjit(rb) + epsi;
                prd  = (gjit >= dm_i);
            else
                % Bernoulli draws using pixel-wise pw
                prd  = rand(numel(rb),1) < pw(rb);
            end

            mB = compute_metrics(prd, y_b);
            BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
            BR.G(b)   = sqrt(max(mB.REC,0)*max(mB.SPEC,0));

            mask_wet_ISPB_b = mask_wet_ISPB_valid(rb);
            mask_dry_OSPB_b = mask_dry_OSPB_valid(rb);

            denom_wet_b = sum( y_b & mask_wet_ISPB_b );
            if denom_wet_b > 0
                TP_wet_b = sum( prd & y_b & mask_wet_ISPB_b );
                rec_wet_b = TP_wet_b / max(denom_wet_b, eps);
            else
                rec_wet_b = NaN;
            end

            denom_dry_b = sum( ~y_b & mask_dry_OSPB_b );
            if denom_dry_b > 0
                TN_dry_b = sum( ~prd & ~y_b & mask_dry_OSPB_b );
                rec_dry_b = TN_dry_b / max(denom_dry_b, eps);
            else
                rec_dry_b = NaN;
            end

            BR.REC_ISPB_wet(b) = rec_wet_b;
            BR.REC_OSPB_dry(b) = rec_dry_b;

            if isfinite(rec_wet_b) && isfinite(rec_dry_b)
                BR.G_ISPBwet_OSPBdry(b) = sqrt(max(rec_wet_b,0) * max(rec_dry_b,0));
            else
                BR.G_ISPBwet_OSPBdry(b) = NaN;
            end

        end

    else
        % -------- Component bootstrap with replacement ----------
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
                rb  = randi(Nf, Nf, 1);
                y_b = yr(rb);

                if use_noisy
                    if isempty(s_all), s_all = sqrt(max(sigM_M.^2 + sigG_M.^2, eps)); end
                    epsi = randn(numel(rb),1) .* (cfg.uncertainty.noise_K .* s_all(rb));
                    gjit = gdif_M_all(V); gjit = gjit(rb) + epsi;
                    prd  = (gjit >= dm_i);
                else
                    prd  = rand(numel(rb),1) < pw(rb);
                end

                mB = compute_metrics(prd, y_b);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0) * max(mB.SPEC,0));

                mask_wet_ISPB_b = mask_wet_ISPB_valid(rb);
                mask_dry_OSPB_b = mask_dry_OSPB_valid(rb);

                denom_wet_b = sum( y_b & mask_wet_ISPB_b );
                if denom_wet_b > 0
                    TP_wet_b = sum( prd & y_b & mask_wet_ISPB_b );
                    rec_wet_b = TP_wet_b / max(denom_wet_b, eps);
                else
                    rec_wet_b = NaN;
                end

                denom_dry_b = sum( ~y_b & mask_dry_OSPB_b );
                if denom_dry_b > 0
                    TN_dry_b = sum( ~prd & ~y_b & mask_dry_OSPB_b );
                    rec_dry_b = TN_dry_b / max(denom_dry_b, eps);
                else
                    rec_dry_b = NaN;
                end

                BR.REC_ISPB_wet(b) = rec_wet_b;
                BR.REC_OSPB_dry(b) = rec_dry_b;

                if isfinite(rec_wet_b) && isfinite(rec_dry_b)
                    BR.G_ISPBwet_OSPBdry(b) = sqrt(max(rec_wet_b,0) * max(rec_dry_b,0));
                else
                    BR.G_ISPBwet_OSPBdry(b) = NaN;
                end
            end
        else
            mfrac = max(0, min(1, cfg.uncertainty.mfrac));
            mapM2V = zeros(numel(Vmask),1); mapM2V(V) = 1:nnz(V);

            % PPS weights
            if cfg.uncertainty.pps_components
                w = cellfun(@numel, GV); sW = sum(w);
                if sW>0, w = w / sW; else, w = []; end
            else
                w = [];
            end

            for b = 1:B
                mC   = max(1, round(mfrac * nC)); % bagging subsample
                if isempty(w)
                    rbC  = randi(nC, mC, 1);            % uniform with replacement
                else
                    % multinomial sampling by CDF
                    u = rand(mC,1); cdf = cumsum(w(:));
                    rbC = arrayfun(@(x) find(cdf>=x,1,'first'), u);
                end

                pick = vertcat(GV{rbC});
                idxV = mapM2V(pick);
                pv   = pw(idxV);
                yv   = yr(idxV);

                if use_noisy
                    if isempty(s_all), s_all = sqrt(max(sigM_M.^2 + sigG_M.^2, eps)); end
                    epsi = randn(numel(idxV),1) .* (cfg.uncertainty.noise_K .* s_all(idxV));
                    gjit = gdif_M_all(V); gjit = gjit(idxV) + epsi;
                    prd  = (gjit >= dm_i);
                else
                    prd  = rand(numel(pv),1) < pv;
                end

                mB = compute_metrics(prd, yv);
                BR.SPEC(b)=mB.SPEC; BR.REC(b)=mB.REC; BR.PREC(b)=mB.PREC; BR.ACC(b)=mB.ACC; BR.F1(b)=mB.F1;
                BR.G(b)   = sqrt(max(mB.REC,0) * max(mB.SPEC,0));

                mask_wet_ISPB_b = mask_wet_ISPB_valid(idxV);
                mask_dry_OSPB_b = mask_dry_OSPB_valid(idxV);

                denom_wet_b = sum( yv & mask_wet_ISPB_b );
                if denom_wet_b > 0
                    TP_wet_b = sum( prd & yv & mask_wet_ISPB_b );
                    rec_wet_b = TP_wet_b / max(denom_wet_b, eps);
                else
                    rec_wet_b = NaN;
                end

                denom_dry_b = sum( ~yv & mask_dry_OSPB_b );
                if denom_dry_b > 0
                    TN_dry_b = sum( ~prd & ~yv & mask_dry_OSPB_b );
                    rec_dry_b = TN_dry_b / max(denom_dry_b, eps);
                else
                    rec_dry_b = NaN;
                end

                BR.REC_ISPB_wet(b) = rec_wet_b;
                BR.REC_OSPB_dry(b) = rec_dry_b;

                if isfinite(rec_wet_b) && isfinite(rec_dry_b)
                    BR.G_ISPBwet_OSPBdry(b) = sqrt(max(rec_wet_b,0) * max(rec_dry_b,0));
                else
                    BR.G_ISPBwet_OSPBdry(b) = NaN;
                end
            end
        end
    end

    % CI computation for base stats + global G
    for s = stat_list
        vR = BR.(s{1});
        CI_raw.(s{1})(i,:)  = prctile(vR,  100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
        CIi_raw.(s{1})(i,:) = prctile(vR,  100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    end
    CI_raw.G(i,:)  = prctile(BR.G,  100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CIi_raw.G(i,:) = prctile(BR.G,  100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    EXP_raw.G(i,1) = sqrt(max(EXP_raw.REC(i),0) * max(EXP_raw.SPEC(i),0));

    % CI computation for regional metrics
    CI_raw.REC_ISPB_wet(i,:)      = prctile(BR.REC_ISPB_wet,      100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CI_raw.REC_OSPB_dry(i,:)      = prctile(BR.REC_OSPB_dry,      100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);
    CI_raw.G_ISPBwet_OSPBdry(i,:) = prctile(BR.G_ISPBwet_OSPBdry, 100*[(1-cfg.uncertainty.band_main)/2, 1 - (1-cfg.uncertainty.band_main)/2]);

    CIi_raw.REC_ISPB_wet(i,:)      = prctile(BR.REC_ISPB_wet,      100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    CIi_raw.REC_OSPB_dry(i,:)      = prctile(BR.REC_OSPB_dry,      100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);
    CIi_raw.G_ISPBwet_OSPBdry(i,:) = prctile(BR.G_ISPBwet_OSPBdry, 100*[(1-cfg.uncertainty.band_inner)/2, 1 - (1-cfg.uncertainty.band_inner)/2]);

end

S.Gstats = struct();
for s = stat_list
    S.Gstats.([s{1} '_raw']) = struct('EXP', EXP_raw.(s{1}), ...
                                      'CI',  CI_raw.(s{1}), ...
                                      'CIi', CIi_raw.(s{1}));
end
S.Gstats.G_raw = struct('EXP', EXP_raw.G, 'CI', CI_raw.G, 'CIi', CIi_raw.G);

S.Gstats.REC_ISPB_wet_raw = struct('EXP', EXP_raw.REC_ISPB_wet, ...
                                   'CI',  CI_raw.REC_ISPB_wet, ...
                                   'CIi', CIi_raw.REC_ISPB_wet);

S.Gstats.REC_OSPB_dry_raw = struct('EXP', EXP_raw.REC_OSPB_dry, ...
                                   'CI',  CI_raw.REC_OSPB_dry, ...
                                   'CIi', CIi_raw.REC_OSPB_dry);

S.Gstats.G_ISPBwet_OSPBdry_raw = struct('EXP', EXP_raw.G_ISPBwet_OSPBdry, ...
                                        'CI',  CI_raw.G_ISPBwet_OSPBdry, ...
                                        'CIi', CIi_raw.G_ISPBwet_OSPBdry);

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

%% σ(GΔ) summaries & CI CSV
S.results_table.sigmaGdel_mean   = S.sigma_gdif_mean;
S.results_table.sigmaGdel_median = S.sigma_gdif_median;
timestamp_sigma = char(datetime('now','Format','yyyyMMdd_HHmmss'));
out_csv_sigma = fullfile(cfg.outdir, sprintf('sigmaGdel_summary%s_%s.csv', region_suffix(cfg), timestamp_sigma));
writetable(S.results_table(:, {'Model','sigmaGdel_mean','sigmaGdel_median'}), out_csv_sigma);
fprintf('Saved sigma(GΔ) summary: %s\n', out_csv_sigma);

okCI = isfield(S,'Gstats') && ...
       isfield(S.Gstats,'SPEC_raw') && isfield(S.Gstats.SPEC_raw,'CI') && ...
       ~isempty(S.Gstats.SPEC_raw.CI) && size(S.Gstats.SPEC_raw.CI,2)==2;

if okCI
  try
    Name = string(S.names(:));
    if isfield(S,'titles') && numel(S.titles)==numel(S.names)
        Title = string(S.titles(:));
    else
        Title = Name;
    end

    % helpers
    getCol = @(fld, sub) i_getCol(S.Gstats, fld, sub, numel(Name));
    getCI  = @(fld, which) i_getCI(S.Gstats, fld, which, numel(Name)); % which = 1(lo) or 2(hi)

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

    haveG = isfield(S.Gstats,'G_raw') && isfield(S.Gstats.G_raw,'EXP') && ...
            ~isempty(S.Gstats.G_raw.EXP) && isfield(S.Gstats.G_raw,'CI') && ...
            ~isempty(S.Gstats.G_raw.CI) && size(S.Gstats.G_raw.CI,2)==2;

    if haveG
      G    = i_pad(S.Gstats.G_raw.EXP(:),  numel(Name));
      G_lo = i_pad(S.Gstats.G_raw.CI(:,1), numel(Name));
      G_hi = i_pad(S.Gstats.G_raw.CI(:,2), numel(Name));
      Tsum.G_raw    = G;
      Tsum.G_raw_lo = G_lo;
      Tsum.G_raw_hi = G_hi;
    end

    % CI summary columns for regional metric
    haveReg = isfield(S.Gstats,'REC_ISPB_wet_raw') && ...
              isfield(S.Gstats.REC_ISPB_wet_raw,'EXP') && ...
              isfield(S.Gstats.REC_ISPB_wet_raw,'CI') && ...
              isfield(S.Gstats,'REC_OSPB_dry_raw') && ...
              isfield(S.Gstats.REC_OSPB_dry_raw,'EXP') && ...
              isfield(S.Gstats.REC_OSPB_dry_raw,'CI') && ...
              isfield(S.Gstats,'G_ISPBwet_OSPBdry_raw') && ...
              isfield(S.Gstats.G_ISPBwet_OSPBdry_raw,'EXP') && ...
              isfield(S.Gstats.G_ISPBwet_OSPBdry_raw,'CI');

    if haveReg
        REC_ISPB_wet    = i_pad(S.Gstats.REC_ISPB_wet_raw.EXP(:), numel(Name));
        REC_ISPB_wet_lo = i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,1), numel(Name));
        REC_ISPB_wet_hi = i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,2), numel(Name));

        REC_OSPB_dry    = i_pad(S.Gstats.REC_OSPB_dry_raw.EXP(:), numel(Name));
        REC_OSPB_dry_lo = i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,1), numel(Name));
        REC_OSPB_dry_hi = i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,2), numel(Name));

        G_reg    = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.EXP(:), numel(Name));
        G_reg_lo = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,1), numel(Name));
        G_reg_hi = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,2), numel(Name));

        Tsum.REC_ISPB_wet_raw    = REC_ISPB_wet;
        Tsum.REC_ISPB_wet_raw_lo = REC_ISPB_wet_lo;
        Tsum.REC_ISPB_wet_raw_hi = REC_ISPB_wet_hi;

        Tsum.REC_OSPB_dry_raw    = REC_OSPB_dry;
        Tsum.REC_OSPB_dry_raw_lo = REC_OSPB_dry_lo;
        Tsum.REC_OSPB_dry_raw_hi = REC_OSPB_dry_hi;

        Tsum.G_ISPBwet_OSPBdry_raw    = G_reg;
        Tsum.G_ISPBwet_OSPBdry_raw_lo = G_reg_lo;
        Tsum.G_ISPBwet_OSPBdry_raw_hi = G_reg_hi;
    end

    timestamp_ci = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    csv_path = fullfile(cfg.outdir, sprintf('ci_summary_%s_%s.csv', region_suffix(cfg), timestamp_ci));
    writetable(Tsum, csv_path);
    fprintf('Saved CI summary CSV: %s\n', csv_path);

  catch ME
    warning('CI/diagnostic readout failed');
  end
end

%% ---- Save lightweight plotting bundle ----    
cfg.artifact_dir    = fullfile(cfg.outdir, 'artifacts');
cfg.ds_plot         = 8;

if ~exist(cfg.artifact_dir,'dir'), mkdir(cfg.artifact_dir); end
ds = max(1, round(cfg.ds_plot)); r = 1:ds:size(S.H,1); c = 1:ds:size(S.H,2);

S_plot = struct();
S_plot.includes = {'X_km','Y_km','H','viz_mask','sink_mask_comp','comp_id', ...
                'roi_mask','y_raw_full','mask_ISPB','mask_OSPB', ...
                'mask_SPB','rect_mask','ds_stride','names','titles','Gstats'};

S_plot.models_ds = struct();
validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);
for ii = 1:numel(S.names)
    fld = validNames{ii};
    if isfield(S.models, fld)
        S_plot.models_ds.(fld) = single(S.models.(fld)(r,c));
    end
end

S_plot.names  = S.names;
if isfield(S,'titles') && numel(S.titles)==numel(S.names), S_plot.titles = S.titles; else, S_plot.titles = S.names; end
S_plot.run_id = cfg.run_id;

S_plot.Gstats   = S.Gstats;
if isfield(S,'model_cvals'), S_plot.model_cvals = S.model_cvals; end
if isfield(S,'model_cvals_median'),  S_plot.model_cvals_median = S.model_cvals_median; end
if isfield(S,'model_quantiles'),     S_plot.model_quantiles    = S.model_quantiles;    end
if isfield(S,'sigma_gdif_mean'),   S_plot.sigma_gdif_mean   = S.sigma_gdif_mean;   end
if isfield(S,'sigma_gdif_median'), S_plot.sigma_gdif_median = S.sigma_gdif_median; end
if isfield(S,'suggested_margin_mean'),   S_plot.suggested_margin_mean   = S.suggested_margin_mean;   end
if isfield(S,'suggested_margin_median'), S_plot.suggested_margin_median = S.suggested_margin_median; end


getCol_plot = @(fld, sub, n) i_getCol(S.Gstats, fld, sub, n);
getCI_plot  = @(fld, which, n) i_getCI(S.Gstats, fld, which, n); % which = 1(lo) | 2(hi)

nM = numel(S.names);

stat_list = {'SPEC','REC','PREC','ACC','F1'};
S_plot.stats = struct();

for k = 1:numel(stat_list)
    fld = [stat_list{k} '_raw'];
    EXP = getCol_plot(fld, 'EXP', nM);
    CI  = [getCI_plot(fld, 1, nM), getCI_plot(fld, 2, nM)];
    CIi = NaN(nM,2);
    if isfield(S.Gstats, fld) && isfield(S.Gstats.(fld), 'CIi') && ~isempty(S.Gstats.(fld).CIi)
        CIi = [i_pad(S.Gstats.(fld).CIi(:,1), nM), i_pad(S.Gstats.(fld).CIi(:,2), nM)];
    end

    S_plot.stats.(stat_list{k}) = struct( ...
        'EXP', EXP, ...
        'CI',  CI, ...
        'CIi', CIi, ...
        'lo',  CI(:,1), ...
        'hi',  CI(:,2), ...
        'lo50', CIi(:,1), ...
        'hi50', CIi(:,2));
end

% Derived G (geometric mean of REC & SPEC): EXP, CI, CIi
S_plot.stats.G = struct('EXP', NaN(nM,1), 'CI', NaN(nM,2), 'CIi', NaN(nM,2), ...
                        'lo', NaN(nM,1), 'hi', NaN(nM,1), 'lo50', NaN(nM,1), 'hi50', NaN(nM,1));
if isfield(S.Gstats,'G_raw')
    S_plot.stats.G.EXP = i_pad(S.Gstats.G_raw.EXP(:), nM);
    if isfield(S.Gstats.G_raw,'CI') && ~isempty(S.Gstats.G_raw.CI)
        S_plot.stats.G.CI  = [i_pad(S.Gstats.G_raw.CI(:,1), nM),  i_pad(S.Gstats.G_raw.CI(:,2), nM)];
        S_plot.stats.G.lo  = S_plot.stats.G.CI(:,1);
        S_plot.stats.G.hi  = S_plot.stats.G.CI(:,2);
    end
    if isfield(S.Gstats.G_raw,'CIi') && ~isempty(S.Gstats.G_raw.CIi)
        S_plot.stats.G.CIi  = [i_pad(S.Gstats.G_raw.CIi(:,1), nM), i_pad(S.Gstats.G_raw.CIi(:,2), nM)];
        S_plot.stats.G.lo50 = S_plot.stats.G.CIi(:,1);
        S_plot.stats.G.hi50 = S_plot.stats.G.CIi(:,2);
    end
else
    % Fallback: compute G EXP from REC/SPEC EXP if needed
    if isfield(S_plot,'stats') && isfield(S_plot.stats,'REC') && isfield(S_plot.stats,'SPEC')
        recE = max(S_plot.stats.REC.EXP(:),0);
        specE = max(S_plot.stats.SPEC.EXP(:),0);
        S_plot.stats.G.EXP = sqrt(recE .* specE);
    end
end

S_plot.stats.REC_ISPB_wet = struct('EXP', NaN(nM,1), 'CI', NaN(nM,2), 'CIi', NaN(nM,2), ...
                                   'lo', NaN(nM,1), 'hi', NaN(nM,1), 'lo50', NaN(nM,1), 'hi50', NaN(nM,1));
S_plot.stats.REC_OSPB_dry = struct('EXP', NaN(nM,1), 'CI', NaN(nM,2), 'CIi', NaN(nM,2), ...
                                   'lo', NaN(nM,1), 'hi', NaN(nM,1), 'lo50', NaN(nM,1), 'hi50', NaN(nM,1));
S_plot.stats.G_ISPBwet_OSPBdry = struct('EXP', NaN(nM,1), 'CI', NaN(nM,2), 'CIi', NaN(nM,2), ...
                                        'lo', NaN(nM,1), 'hi', NaN(nM,1), 'lo50', NaN(nM,1), 'hi50', NaN(nM,1));

if isfield(S.Gstats,'REC_ISPB_wet_raw')
    S_plot.stats.REC_ISPB_wet.EXP = i_pad(S.Gstats.REC_ISPB_wet_raw.EXP(:), nM);
    if isfield(S.Gstats.REC_ISPB_wet_raw,'CI') && ~isempty(S.Gstats.REC_ISPB_wet_raw.CI)
        S_plot.stats.REC_ISPB_wet.CI  = [i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,1), nM), i_pad(S.Gstats.REC_ISPB_wet_raw.CI(:,2), nM)];
        S_plot.stats.REC_ISPB_wet.lo  = S_plot.stats.REC_ISPB_wet.CI(:,1);
        S_plot.stats.REC_ISPB_wet.hi  = S_plot.stats.REC_ISPB_wet.CI(:,2);
    end
    if isfield(S.Gstats.REC_ISPB_wet_raw,'CIi') && ~isempty(S.Gstats.REC_ISPB_wet_raw.CIi)
        S_plot.stats.REC_ISPB_wet.CIi  = [i_pad(S.Gstats.REC_ISPB_wet_raw.CIi(:,1), nM), i_pad(S.Gstats.REC_ISPB_wet_raw.CIi(:,2), nM)];
        S_plot.stats.REC_ISPB_wet.lo50 = S_plot.stats.REC_ISPB_wet.CIi(:,1);
        S_plot.stats.REC_ISPB_wet.hi50 = S_plot.stats.REC_ISPB_wet.CIi(:,2);
    end
end

if isfield(S.Gstats,'REC_OSPB_dry_raw')
    S_plot.stats.REC_OSPB_dry.EXP = i_pad(S.Gstats.REC_OSPB_dry_raw.EXP(:), nM);
    if isfield(S.Gstats.REC_OSPB_dry_raw,'CI') && ~isempty(S.Gstats.REC_OSPB_dry_raw.CI)
        S_plot.stats.REC_OSPB_dry.CI  = [i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,1), nM), i_pad(S.Gstats.REC_OSPB_dry_raw.CI(:,2), nM)];
        S_plot.stats.REC_OSPB_dry.lo  = S_plot.stats.REC_OSPB_dry.CI(:,1);
        S_plot.stats.REC_OSPB_dry.hi  = S_plot.stats.REC_OSPB_dry.CI(:,2);
    end
    if isfield(S.Gstats.REC_OSPB_dry_raw,'CIi') && ~isempty(S.Gstats.REC_OSPB_dry_raw.CIi)
        S_plot.stats.REC_OSPB_dry.CIi  = [i_pad(S.Gstats.REC_OSPB_dry_raw.CIi(:,1), nM), i_pad(S.Gstats.REC_OSPB_dry_raw.CIi(:,2), nM)];
        S_plot.stats.REC_OSPB_dry.lo50 = S_plot.stats.REC_OSPB_dry.CIi(:,1);
        S_plot.stats.REC_OSPB_dry.hi50 = S_plot.stats.REC_OSPB_dry.CIi(:,2);
    end
end

if isfield(S.Gstats,'G_ISPBwet_OSPBdry_raw')
    S_plot.stats.G_ISPBwet_OSPBdry.EXP = i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.EXP(:), nM);
    if isfield(S.Gstats.G_ISPBwet_OSPBdry_raw,'CI') && ~isempty(S.Gstats.G_ISPBwet_OSPBdry_raw.CI)
        S_plot.stats.G_ISPBwet_OSPBdry.CI  = [i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,1), nM), i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CI(:,2), nM)];
        S_plot.stats.G_ISPBwet_OSPBdry.lo  = S_plot.stats.G_ISPBwet_OSPBdry.CI(:,1);
        S_plot.stats.G_ISPBwet_OSPBdry.hi  = S_plot.stats.G_ISPBwet_OSPBdry.CI(:,2);
    end
    if isfield(S.Gstats.G_ISPBwet_OSPBdry_raw,'CIi') && ~isempty(S.Gstats.G_ISPBwet_OSPBdry_raw.CIi)
        S_plot.stats.G_ISPBwet_OSPBdry.CIi  = [i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CIi(:,1), nM), i_pad(S.Gstats.G_ISPBwet_OSPBdry_raw.CIi(:,2), nM)];
        S_plot.stats.G_ISPBwet_OSPBdry.lo50 = S_plot.stats.G_ISPBwet_OSPBdry.CIi(:,1);
        S_plot.stats.G_ISPBwet_OSPBdry.hi50 = S_plot.stats.G_ISPBwet_OSPBdry.CIi(:,2);
    end
end

S_plot.X_km     = single(S.X_km(r,c));
S_plot.Y_km     = single(S.Y_km(r,c));
S_plot.H        = single(S.H(r,c));
S_plot.viz_mask = logical(S.viz_mask(r,c));

S_plot.sink_mask_comp = logical(S.sink_mask_comp(r,c));
S_plot.comp_id        = uint32(S.comp_id(r,c));
S_plot.roi_mask       = logical(S.roi_mask(r,c));
S_plot.y_raw_full     = logical(S.y_raw_full(r,c));

if isfield(S,'mask_ISPB'), S_plot.mask_ISPB = logical(S.mask_ISPB(r,c)); end
if isfield(S,'mask_OSPB'), S_plot.mask_OSPB = logical(S.mask_OSPB(r,c)); end
if isfield(S,'mask_SPB'),  S_plot.mask_SPB  = logical(S.mask_SPB(r,c));  end
if isfield(S,'rect_mask'), S_plot.rect_mask = logical(S.rect_mask(r,c)); end

S_plot.ds_stride = uint16(ds);
if isfield(cfg,'decision_margin_mode'), S_plot.margin_used_mode = string(cfg.decision_margin_mode); else, S_plot.margin_used_mode = "fixed"; end

cfg_plot = struct();
cfg_plot.outdir    = cfg.outdir;
cfg_plot.run_id    = cfg.run_id;
cfg_plot.font_size = cfg.font_size;
cfg_plot.overlay_alpha = cfg.overlay_alpha;
if isfield(cfg,'region'), cfg_plot.region = cfg.region; else, cfg_plot.region = struct('mode','ALL'); end
cfg_plot.flat = struct('interp_enable',true,'interp_points',50,'interp_method','spline');

stamp  = char(datetime('now','Format','yyyyMMdd_HHmmss'));
f_plot = fullfile(cfg.artifact_dir, sprintf('S_plot_small_%s.mat', stamp));
save(f_plot, 'S_plot', 'cfg_plot', '-v7');
fprintf('[artifacts] Saved lightweight plot file: %s\n', f_plot);

cache_dir = fullfile(cfg.outdir, 'gdif_cache'); if ~exist(cache_dir,'dir'), mkdir(cache_dir); end
is_flat = startsWith(S.names,'Flat_'); validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);
for ii = 1:numel(S.names)
    if is_flat(ii), continue; end
    name = S.names{ii}; fld = validNames{ii};
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

end 

%% helpers
function s = region_suffix(cfg)
    s = '';
    try
        if isfield(cfg,'region') && isfield(cfg.region,'mode') && ~isempty(cfg.region.mode)
            s = ['_' upper(char(cfg.region.mode))];
        end
    catch
    end
end

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
        if isfield(S,'sigma_gdif_mean'), nf_mask = true(numel(S.sigma_gdif_mean),1);
        else, nf_mask = true(1,1);
        end
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
            if haveMean && exist('iModel','var') && iModel>=1 && iModel<=numel(S.sigma_gdif_mean) && isfinite(S.sigma_gdif_mean(iModel))
                dm = S.sigma_gdif_mean(iModel);
            else, dm = cfg.decision_margin; 
            end
        case "auto_per_model_median"
            if haveMedian && exist('iModel','var') && iModel>=1 && iModel<=numel(S.sigma_gdif_median) && isfinite(S.sigma_gdif_median(iModel))
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
