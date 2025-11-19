function R = make_synthetic_GHF_MCMC_new(cfg)
% MAKE_SYNTHETIC_GHF_MCMC (streamlined, physics-gated) — whole-grid mean + spec-valid median
%   Single-chain MCMC maximizing performance metric:
%       G_ISPBwet_OSPBdry = sqrt( REC(wet in ISPB) * REC(dry in OSPB) )
%
%   Records:
%     - G_mu_grid:  whole-grid mean(G) after gating (finite pixels)
%     - G_median_specvalid: median(G) over all spec-valid pixels (Q ~= 9999), independent of ROI
%
% Outputs struct R with fields:
%   mode, chain (table), keep (table), best (struct), acc_rate, run_id, etc.

fprintf('\n=============================\n');
fprintf('[make_synthetic_GHF_MCMC] Start (%s)\n', datetime('now'));
fprintf('=============================\n');

%% ---------------- Defaults ----------------
if nargin<1, cfg = struct(); end
cfg = set_default(cfg, 'template_mode',   'gradient');   % 'gradient' | 'randn-smoothed'
cfg = set_default(cfg, 'sample_angle',    true);
cfg = set_default(cfg, 'n_steps',         12000);
cfg = set_default(cfg, 'burn_frac',       0.25);
cfg = set_default(cfg, 'thin',            1);
cfg = set_default(cfg, 'beta',            6);            % posterior scale for regional G-mean
cfg = set_default(cfg, 'step_mu',         30);
cfg = set_default(cfg, 'step_het',        30);
cfg = set_default(cfg, 'step_theta',      30);
cfg = set_default(cfg, 'mu_range',        [20 150]);
cfg = set_default(cfg, 'het_range',       [-150 150]);   % A (gradient) or sigma (noise)
cfg = set_default(cfg, 'theta_range_deg', [-90 90]);
cfg = set_default(cfg, 'prior_type',      'softbox');    % 'softbox' | 'gaussian'
cfg = set_default(cfg, 'soft_k',          3);
cfg = set_default(cfg, 'seed',            42);
cfg = set_default(cfg, 'outdir',          'figs_out');
cfg = set_default(cfg, 'region_mode',     'ALL');        % recommended: 'ALL' so ISPB+OSPB both present
cfg = set_default(cfg, 'verbose',         true);
cfg = set_default(cfg, 'save_maps',       true);
cfg = set_default(cfg, 'map_every',       1000);
cfg = set_default(cfg, 'map_downsample',  1);
cfg = set_default(cfg, 'map_clim',        struct('mode','symmetric','k',2,'pad',0));
cfg = set_default(cfg, 'grad_angle_deg0', 30);           % if sample_angle=false
cfg = set_default(cfg, 'roi_x_range',     [0 700]*1000);
cfg = set_default(cfg, 'roi_y_range',     [-200 350]*1000);
cfg = set_default(cfg, 'proxy_gate',      0);            % skip full eval if proxy metric below this
cfg = set_default(cfg, 'cache_round',     [0.25, 0.5, 2]); % rounding for [mu, A/sigma, theta]

% --- Physics gating knobs (unified names) ---
cfg = set_default(cfg, 'floor_G',         20);     % lower bound (mW m^-2)
cfg = set_default(cfg, 'cap_G',           100);    % upper bound (mW m^-2), [] => none
cfg = set_default(cfg, 'enforce_mean',    true);   % recenter to target mu after floor/cap
cfg = set_default(cfg, 'A_nonneg',        false);  % restrict A >= 0 (absorb sign into theta)
cfg = set_default(cfg, 'penalize_clip',   true);   % soft penalty for expected clipping frac
cfg = set_default(cfg, 'clip_k',          0.15);   % softness for clipping penalty
cfg = set_default(cfg, 'skip_bootstrap',  true);   % avoid heavy side-effects during eval

% Sanity checks
assert(numel(cfg.mu_range)==2 && cfg.mu_range(1)<cfg.mu_range(2), 'mu_range invalid');
assert(numel(cfg.het_range)==2 && cfg.het_range(1)<cfg.het_range(2), 'het_range invalid');
if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
    assert(numel(cfg.theta_range_deg)==2 && cfg.theta_range_deg(1)<cfg.theta_range_deg(2), 'theta_range invalid');
end

% Output dirs / RNG
run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
artdir = fullfile(cfg.outdir,'artifacts'); if ~exist(artdir,'dir'), mkdir(artdir); end
mapdir = fullfile(cfg.outdir,'maps');      if ~exist(mapdir,'dir'), mkdir(mapdir); end
rng(cfg.seed);

if cfg.verbose
    fprintf('[config] template=%s | steps=%d | burn=%.0f%%%% | thin=%d | beta=%.2f | proxy_gate=%.3f\n', ...
        cfg.template_mode, cfg.n_steps, 100*cfg.burn_frac, cfg.thin, cfg.beta, cfg.proxy_gate);
end

%% ---------------- Grid & static state ----------------
if isfield(cfg,'S_static') && ~isempty(cfg.S_static)
    S_static = cfg.S_static;  % allow caller to inject a prebuilt light state
    assert(isfield(S_static,'Xgrid') && isfield(S_static,'Ygrid'), 'S_static missing Xgrid/Ygrid');
else
    D = load('datasets_for_gmin/coldex_icethk.mat');
    assert(isfield(D,'Xgrid') && isfield(D,'Ygrid'), 'Xgrid/Ygrid missing in coldex_icethk.mat');
    S_static = load('datasets_for_gmin/gmin_data.mat','S');
    S_static = S_static.S;   % contains Xgrid,Ygrid,Q,H,icevel,[Gmin]
end
X = single(S_static.Xgrid); Y = single(S_static.Ygrid);

% Attach Gmin if stored separately
if ~isfield(S_static,'Gmin')
    try L2 = load('datasets_for_gmin/gmin_data.mat','Gmin'); if isfield(L2,'Gmin'), S_static.Gmin = L2.Gmin; end; catch, end
end

% Region masks and ROI
Ms = load('datasets_for_gmin/spb_masks.mat');
S_static.mask_ISPB = false(size(X)); S_static.mask_OSPB = S_static.mask_ISPB; S_static.mask_SPB = S_static.mask_ISPB;
S_static.mask_ISPB(Ms.rmin:Ms.rmax, Ms.cmin:Ms.cmax) = logical(full(Ms.mask_ISPB_c));
S_static.mask_OSPB(Ms.rmin:Ms.rmax, Ms.cmin:Ms.cmax) = logical(full(Ms.mask_OSPB_c));
S_static.mask_SPB( Ms.rmin:Ms.rmax, Ms.cmin:Ms.cmax) = logical(full(Ms.mask_SPB_c));
S_static.rect_mask = logical(Ms.rect_mask);

spec_thr = 0.2; v_keep = 10;
if ~isfield(S_static,'spec_invalid') && isfield(S_static,'Q')
    S_static.spec_invalid = (S_static.Q == 9999);
end
S_static.spec_valid_mask = isfinite(S_static.Q) & ~S_static.spec_invalid; % persistent spec-valid

valid_mask = isfinite(S_static.Q) & ~S_static.spec_invalid & isfinite(S_static.H) & isfinite(S_static.icevel);

switch lower(cfg.region_mode)
  case 'ispb',  REG_MASK = S_static.mask_ISPB & S_static.rect_mask;
  case 'ospb',  REG_MASK = S_static.mask_OSPB & S_static.rect_mask;
  case 'spb',   REG_MASK = S_static.mask_SPB  & S_static.rect_mask;
  case 'all',   REG_MASK = S_static.rect_mask;
  otherwise, error('region_mode invalid.');
end
slow_mask         = valid_mask & ~(S_static.icevel > v_keep);
S_static.roi_mask = slow_mask & REG_MASK;

% Optional coordinate clamp
if ~isempty(cfg.roi_x_range) || ~isempty(cfg.roi_y_range)
    coord_mask = true(size(X));
    if ~isempty(cfg.roi_x_range), coord_mask = coord_mask & X >= cfg.roi_x_range(1) & X <= cfg.roi_x_range(2); end
    if ~isempty(cfg.roi_y_range), coord_mask = coord_mask & Y >= cfg.roi_y_range(1) & Y <= cfg.roi_y_range(2); end
    S_static.roi_mask = S_static.roi_mask & coord_mask;
end

% Freeze ROI vectors (no allocations inside loop)
idx = find(S_static.roi_mask);
S_static.roi_idx   = idx;
spec_roi           = S_static.Q(idx);
S_static.y_raw_vec = isfinite(spec_roi) & (spec_roi ~= 9999) & (spec_roi > spec_thr);

% Region labels on ROI-vector (for new metric)
if isfield(S_static,'mask_ISPB') && isfield(S_static,'mask_OSPB')
    S_static.roi_ISPB = logical(S_static.mask_ISPB(idx));  % same length as roi_idx
    S_static.roi_OSPB = logical(S_static.mask_OSPB(idx));
else
    S_static.roi_ISPB = false(numel(idx),1);
    S_static.roi_OSPB = false(numel(idx),1);
end

% Sanity prints
n_roi = numel(idx); n_pos = nnz(S_static.y_raw_vec); n_neg = n_roi - n_pos;
fprintf('[check] ROI=%d | pos=%d (%.2f%%) | neg=%d (%.2f%%)\n', n_roi, n_pos, 100*n_pos/max(1,n_roi), n_neg, 100*n_neg/max(1,n_roi));
assert(isfield(S_static,'Gmin') && isequal(size(S_static.Gmin), size(X)), 'Gmin missing or wrong size');
roi = S_static.roi_mask;
pctGminFinite = 100*nnz(isfinite(S_static.Gmin(roi)))/max(1,nnz(roi));
fprintf('[check] Gmin finite in ROI: %.1f%%%%\n', pctGminFinite);

%% ---------------- Precompute ramps (gradient mode) ----------------
rampCache = containers.Map('KeyType','char','ValueType','any');
if strcmpi(cfg.template_mode,'gradient')
    th_lo = ceil(cfg.theta_range_deg(1)); th_hi = floor(cfg.theta_range_deg(2));
    for th = th_lo:th_hi
        rampCache(sprintf('th_%d', th)) = unit_gradient_angle(X, Y, th); % single
    end
end

%% ---------------- Initial params & quick probes ----------------
mu0  = 40;
het0 = mean(cfg.het_range);
params0 = [mu0, het0];
if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
    params0 = [mu0, het0, 0]; % center angle
end

probe_list = [mu0, het0, 0; mu0, max(cfg.het_range)*0.8, 0; mu0, min(cfg.het_range)*0.2, 0];
probe_list = probe_list(:,1:numel(params0));
for i=1:size(probe_list,1)
    [Gtmp, pst] = build_G_from_params(X, Y, cfg, probe_list(i,:), rampCache);
    [g_reg, r_IS, r_OS, mu_grid, med_spec]  = evaluate_params(Gtmp, S_static, cfg, true);
    fprintf('[probe %d] %s -> G_ISPBwet×OSPBdry=%.3f (REC_ISPB_wet=%.3f, REC_OSPB_dry=%.3f) | <G>=%.2f | med_spec(G)=%.2f\n', ...
        i, fmtv(probe_list(i,:)), g_reg, r_IS, r_OS, mu_grid, med_spec);
    print_phys_stats('probe', pst);
end

[G0, pst0] = build_G_from_params(X, Y, cfg, params0, rampCache);
[gmean0, rec0, spec0, mu_grid0, med_spec0] = evaluate_params(G0, S_static, cfg, true);
lp0   = log_prior(params0, cfg);
post0 = cfg.beta * clamp01(gmean0) + lp0;

% Storage
K = cfg.n_steps; Ddim = numel(params0); params = params0;
chain_params = nan(K, Ddim); chain_score = nan(K, 1);
chain_lp = nan(K, 1); chain_post = nan(K, 1); chain_acc = false(K,1);
chain_gmu_grid   = nan(K,1);   % whole-grid mean(G) after gating
chain_gmed_spec  = nan(K,1);   % spec-valid median(G) after gating

best = struct('params',params,'gmean',gmean0,'rec',rec0,'spec',spec0,'post',post0, ...
              'iter',1,'mu_grid',mu_grid0,'med_spec',med_spec0);

fprintf(['[mcmc] Start @ params=%s | G_ISPBwet×OSPBdry=%.3f ' ...
         '(REC_ISPB_wet=%.3f, REC_OSPB_dry=%.3f) | post=%.3f | <G>=%.2f | med_spec(G)=%.2f\n'], ...
        fmtv(params), gmean0, rec0, spec0, post0, mu_grid0, med_spec0);
print_phys_stats('start', pst0);

% Memo cache (full eval only)
% Cache layout: [g_reg, REC_ISPB_wet, REC_OSPB_dry, mu_grid, med_spec]
ScoreCache = containers.Map('KeyType','char','ValueType','any');
accepts = 0; last_map_iter = 0;

%% ---------------- MCMC loop ----------------
for k = 1:K
    % Propose (random-walk)
    prop = params;
    prop(1) = params(1) + cfg.step_mu  * randn;      % mu
    prop(2) = params(2) + cfg.step_het * randn;      % A or sigma
    if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
        prop(3) = wrap_theta(params(3) + cfg.step_theta * randn, cfg.theta_range_deg);
    end

    key = cache_key(prop, cfg.cache_round, cfg);
    used_cache = false;

    % --- Delayed-acceptance proxy using regional metric ---
    do_full_eval = true;
    if cfg.proxy_gate > 0
        [Gproxy, ~] = build_G_from_params(X, Y, cfg, prop, rampCache);
        p_gm = proxy_gmean(Gproxy, S_static);  % now returns regional G-mean
        if p_gm < cfg.proxy_gate
            do_full_eval = false; % early reject unless prior boosts it (rare)
        end
    end

    if isKey(ScoreCache, key)
        ss = ScoreCache(key);
        gmean1    = ss(1);
        rec1      = ss(2);
        spec1     = ss(3);
        mu_grid1  = ss(4);
        med_spec1 = ss(5);
        used_cache = true; pst1 = [];
    elseif do_full_eval
        [G1, pst1] = build_G_from_params(X, Y, cfg, prop, rampCache);
        [gmean1, rec1, spec1, mu_grid1, med_spec1] = evaluate_params(G1, S_static, cfg, false);
        ScoreCache(key) = [gmean1, rec1, spec1, mu_grid1, med_spec1];
    else
        gmean1 = 0; rec1=0; spec1=1; mu_grid1 = NaN; med_spec1 = NaN; pst1=[];
    end

    lp1   = log_prior(prop, cfg);
    post1 = cfg.beta * clamp01(gmean1) + lp1;

    % MH accept
    dpost = post1 - post0;

    if (dpost >= 0) || (rand < exp(dpost))
        params = prop; gmean0 = gmean1; rec0 = rec1; spec0 = spec1; lp0 = lp1; post0 = post1;
        mu_grid0  = mu_grid1;
        med_spec0 = med_spec1;   % advance tracked stats with state
        chain_acc(k) = true; accepts = accepts + 1;

        if isempty(pst1) && do_full_eval
            [~, pst1] = build_G_from_params(X, Y, cfg, params, rampCache);
        end
        if ~isempty(pst1), print_phys_stats(sprintf('acc %4d',k), pst1); end

        if gmean0 > best.gmean
            best = struct('params',params,'gmean',gmean0,'rec',rec0,'spec',spec0,'post',post0, ...
                          'iter',k,'mu_grid',mu_grid0,'med_spec',med_spec0);
            if cfg.save_maps && (k - last_map_iter >= cfg.map_every || k==1)
                [GbestIter, ~] = build_G_from_params(X, Y, cfg, params, rampCache);
                save_one_map(GbestIter, X, Y, cfg, mapdir, sprintf('BEST_iter%d',k)); % keep closed during run
                last_map_iter = k;
            end
        end
    end

    chain_params(k,:) = params;
    chain_score(k)    = gmean0;      % regional G-mean
    chain_lp(k)       = lp0;
    chain_post(k)     = post0;
    chain_gmu_grid(k)  = mu_grid0;
    chain_gmed_spec(k) = med_spec0;

    if cfg.verbose && mod(k, max(1,round(K/20)))==0
        fprintf(['[mcmc] %4d/%4d | acc=%.1f%%%% | G_ISPBwet×OSPBdry=%.3f ' ...
                 '| post=%.3f | <G>=%.2f | med_spec(G)=%.2f | params=%s%s\n'], ...
            k, K, 100*accepts/k, gmean0, post0, mu_grid0, med_spec0, fmtv(params), ternary(used_cache,' [cache]',''));
    end
end

% --- Final best map (kept open) ---
[GbestFinal, ~] = build_G_from_params(X, Y, cfg, best.params, rampCache);
save_one_map(GbestFinal, X, Y, cfg, mapdir, sprintf('BEST_final_iter%d', best.iter), [], true); % keep_open = true

% Also emit a GeoTIFF for QGIS
qgis_save_geotiff(GbestFinal, X, Y, fullfile(mapdir, sprintf('ghf_best_%s.tif', run_id)), -9999);

acc_rate = accepts / K;
fprintf(['[mcmc] Done. Acceptance rate = %.1f%%%% | Best G_ISPBwet×OSPBdry=%.3f @ %s ' ...
         '(REC_ISPB_wet=%.3f, REC_OSPB_dry=%.3f) | <G>=%.2f | med_spec(G)=%.2f\n'], ...
    100*acc_rate, best.gmean, fmtv(best.params), best.rec, best.spec, best.mu_grid, best.med_spec);

%% ---------------- Post & save ----------------
burn = max(0, round(cfg.burn_frac * K));
keep_idx    = (burn+1):cfg.thin:K;
keep_params = chain_params(keep_idx,:);
keep_score  = chain_score(keep_idx);
keep_post   = chain_post(keep_idx);
keep_mu     = chain_gmu_grid(keep_idx);
keep_medS   = chain_gmed_spec(keep_idx);

names = param_names(cfg);

% NOTE: "Gmean" now means regional G_ISPBwet×OSPBdry
T_all  = array2table(chain_params, 'VariableNames', names);
T_all.Gmean               = chain_score;
T_all.LogPrior            = chain_lp;
T_all.Posterior           = chain_post;
T_all.Accepted            = chain_acc;
T_all.G_mu_grid           = chain_gmu_grid;     % mean over all finite pixels
T_all.G_median_specvalid  = chain_gmed_spec;    % median over spec-valid pixels (Q ~= 9999)

T_keep = array2table(keep_params, 'VariableNames', names);
T_keep.Gmean               = keep_score;
T_keep.Posterior           = keep_post;
T_keep.G_mu_grid           = keep_mu;
T_keep.G_median_specvalid  = keep_medS;

plot_trace(chain_score, 'G_ISPBwet×OSPBdry', cfg, run_id, cfg.outdir);
for d=1:Ddim, plot_trace(chain_params(:,d), names{d}, cfg, run_id, cfg.outdir); end
plot_posterior(keep_params, names, cfg, run_id, cfg.outdir);
plot_pairs(keep_params, names, cfg, run_id, cfg.outdir);

R = struct();
R.mode          = 'mcmc';
R.chain         = T_all;
R.keep          = T_keep;
R.best          = best;
R.acc_rate      = acc_rate;
R.run_id        = run_id;
R.outdir        = cfg.outdir;
R.region_mode   = cfg.region_mode;
R.seed          = cfg.seed;
R.template_mode = cfg.template_mode;
R.sample_angle  = cfg.sample_angle;
R.mu_range      = cfg.mu_range;
R.het_range     = cfg.het_range;
if strcmpi(cfg.template_mode,'gradient'), R.theta_range_deg = cfg.theta_range_deg; end
R.roi_x_range   = cfg.roi_x_range;
R.roi_y_range   = cfg.roi_y_range;

outfile = fullfile(artdir, sprintf('mcmc_results_%s.mat', run_id));
save(outfile, '-struct', 'R');
fprintf('[save] Results saved: %s\n', outfile);
fprintf('[done] Complete. Artifacts in %s\n', artdir);
fprintf('=============================\n\n');
end

% ============================== helpers ==============================

function cfg = set_default(cfg, f, v)
if ~isfield(cfg,f) || isempty(cfg.(f)), cfg.(f) = v; end
end

function key = cache_key(p, roundv, cfg)
p = p(:).';
% respect A_nonneg in key to align with build-time enforcement
if strcmpi(cfg.template_mode,'gradient')
    Aeff = p(2);
    if isfield(cfg,'A_nonneg') && cfg.A_nonneg, Aeff = max(Aeff, 0); end
    if cfg.sample_angle
        p = [p(1), Aeff, p(3)];
    else
        p = [p(1), Aeff];
    end
end
% round
if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
    p = [round(p(1)/roundv(1))*roundv(1), ...
         round(p(2)/roundv(2))*roundv(2), ...
         round(p(3)/roundv(3))*roundv(3)];
elseif numel(p)>=2
    p = [round(p(1)/roundv(1))*roundv(1), ...
         round(p(2)/roundv(2))*roundv(2)];
end
% format key
if numel(p)==3, key = sprintf('%.2f|%.2f|%.1f', p);
else,           key = sprintf('%.2f|%.2f', p);
end
end

function [G, pst] = build_G_from_params(X, Y, cfg, params, rampCache)
switch lower(cfg.template_mode)
case 'gradient'
    mu = single(params(1));
    A  = single(params(2));
    if isfield(cfg,'A_nonneg') && cfg.A_nonneg, A = max(A, 0); end
    if cfg.sample_angle, th = params(3); else, th = cfg.grad_angle_deg0; end
    key = sprintf('th_%d', round(th));
    if isKey(rampCache, key), U = rampCache(key); else, U = unit_gradient_angle(X, Y, th); end
    G = mu + A * U;  % raw field
case 'randn-smoothed'
    mu = single(params(1)); sig = max(single(params(2)), 0);
    Z  = imgaussfilt(randn(size(X),'single'), 15);
    Z  = standardize_Z(Z);
    G  = mu + sig * Z;
otherwise
    error('Unknown template_mode: %s', cfg.template_mode);
end

% ---- Physical gating ----
Gm_before = mean(G(:), 'omitnan');
if ~isempty(cfg.floor_G), G = max(G, single(cfg.floor_G)); end
if ~isempty(cfg.cap_G),   G = min(G, single(cfg.cap_G));   end
if cfg.enforce_mean
    G = G - Gm_before + single(params(1));
    if ~isempty(cfg.floor_G), G = max(G, single(cfg.floor_G)); end
    if ~isempty(cfg.cap_G),   G = min(G, single(cfg.cap_G));   end
end

% stats (whole-grid)
[f_floor, f_cap] = physics_clip_stats(G, cfg.floor_G, cfg.cap_G, []);
pst = struct('min_after',min(G(:),[],'omitnan'), 'max_after',max(G(:),[],'omitnan'), ...
             'mean_after',mean(G(:),'omitnan'), 'f_floor',f_floor, 'f_cap',f_cap);
end

function print_phys_stats(prefix, pst)
try
    fprintf('[phys|%s] min=%.2f max=%.2f mean=%.2f | floor%%=%.2f cap%%=%.2f\n', ...
        prefix, pst.min_after, pst.max_after, pst.mean_after, 100*pst.f_floor, 100*pst.f_cap);
catch
    % ignore if struct incomplete
end
end

function g = proxy_gmean(G, S_static)
% PROXY_GMEAN
%   Cheap surrogate for the regional metric using a median-threshold
%   classifier on raw G within ROI.
idx = S_static.roi_idx;
if isempty(idx) || ~isfield(S_static,'y_raw_vec') || isempty(S_static.y_raw_vec)
    g = 0; return;
end

y = logical(S_static.y_raw_vec(:));
Gi = double(G(idx));
roi_ISPB = S_static.roi_ISPB(:);
roi_OSPB = S_static.roi_OSPB(:);

bad = ~isfinite(Gi) | ~isfinite(y);
if any(bad)
    Gi(bad)      = [];
    y(bad)       = [];
    roi_ISPB(bad)= [];
    roi_OSPB(bad)= [];
end
if isempty(Gi) || numel(y)~=numel(Gi)
    g = 0; return;
end

thr  = median(Gi,'omitnan');
yhat = Gi >= thr;

% regional recalls
pos_IS = y & roi_ISPB;
neg_OS = ~y & roi_OSPB;

N_wet_ISPB  = nnz(pos_IS);
N_dry_OSPB  = nnz(neg_OS);

if N_wet_ISPB==0 || N_dry_OSPB==0
    % fallback: global confusion
    [rec_all, spec_all] = local_confusion(y, yhat);
    g = sqrt(rec_all*spec_all);
else
    REC_ISPB_wet  = nnz(pos_IS & yhat) / max(1,N_wet_ISPB);
    TN_OSPB       = nnz(neg_OS & ~yhat);
    REC_OSPB_dry  = TN_OSPB / max(1,N_dry_OSPB);
    g             = sqrt(max(REC_ISPB_wet,0)*max(REC_OSPB_dry,0));
end

g = clamp01(g);
end

function [gmean, REC_ISPB_wet, REC_OSPB_dry, mu_grid, med_spec] = evaluate_params(G, S_static, cfg, verbose)
% Robust evaluator using GDelta = G - Gmin inside frozen ROI.
% Returns:
%   gmean        : regional G_ISPBwet×OSPBdry
%   REC_ISPB_wet : recall of wet sinks in ISPB
%   REC_OSPB_dry : recall of dry sinks in OSPB (true-negative recall)
%   mu_grid      : whole-grid mean(G) after gating (finite pixels)
%   med_spec     : median(G) over all spec-valid pixels (Q ~= 9999), independent of ROI

idx = S_static.roi_idx; y = logical(S_static.y_raw_vec(:));
assert(numel(idx)==numel(y), 'ROI/label length mismatch.');

% Whole-grid mean(G) (post-gating)
Gf = G(isfinite(G));
if isempty(Gf)
    mu_grid = NaN;
else
    mu_grid = mean(Gf,'omitnan');
end

% Spec-valid median(G) (post-gating), independent of ROI
if isfield(S_static,'spec_valid_mask') && ~isempty(S_static.spec_valid_mask)
    msv = S_static.spec_valid_mask & isfinite(G);
else
    msv = isfinite(S_static.Q) & (S_static.Q ~= 9999) & isfinite(G);
end
gsv = double(G(msv));
if isempty(gsv)
    med_spec = NaN;
else
    med_spec = median(gsv,'omitnan');
end

% ROI evaluation on G - Gmin
g  = double(G(idx));
gm = double(S_static.Gmin(idx));

bad = ~isfinite(g) | ~isfinite(gm);
if any(bad)
    g(bad)  = [];
    gm(bad) = [];
    y(bad)  = [];
end

n    = numel(y);
npos = nnz(y);
nneg = n - npos;

if n==0 || npos==0 || nneg==0
    % fallback to proxy-like global G-mean
    thr  = median(double(G(idx)),'omitnan');
    Gi   = double(G(idx));
    bad2 = ~isfinite(Gi);
    if any(bad2)
        Gi(bad2) = [];
        y2 = y; y2(bad2) = [];
    else
        y2 = y;
    end
    yhat = Gi >= thr;
    [rec_all, spec_all] = local_confusion(y2, yhat);
    gmean          = sqrt(rec_all*spec_all);
    REC_ISPB_wet   = NaN;
    REC_OSPB_dry   = NaN;
    if verbose
        fprintf('[triage] Degenerate global classes; proxy global G-mean used. REC=%.3f SPEC=%.3f Gmn=%.3f\n', ...
            rec_all, spec_all, gmean);
    end
    return;
end

% decision rule on GDelta
yhat = (g - gm) >= 0;

% global confusion (for diagnostics / fallback)
[rec_all, spec_all] = local_confusion(y, yhat);
G_global = sqrt(rec_all*spec_all);

% regional labels (after same "bad" filter)
roi_ISPB = S_static.roi_ISPB(:);
roi_OSPB = S_static.roi_OSPB(:);
if any(bad)
    roi_ISPB(bad) = [];
    roi_OSPB(bad) = [];
end

pos_IS    = y  & roi_ISPB;
neg_OS    = ~y & roi_OSPB;
N_wet_ISPB = nnz(pos_IS);
N_dry_OSPB = nnz(neg_OS);

if N_wet_ISPB==0 || N_dry_OSPB==0
    % if one region has no relevant class, fall back to global as objective
    REC_ISPB_wet  = NaN;
    REC_OSPB_dry  = NaN;
    gmean         = G_global;
    if verbose
        fprintf(['[eval:Gmin] Regional counts degenerate (Nwet_ISPB=%d, Ndry_OSPB=%d); ' ...
                 'using global G-mean=%.3f (REC_all=%.3f, SPEC_all=%.3f)\n'], ...
                 N_wet_ISPB, N_dry_OSPB, gmean, rec_all, spec_all);
    end
    return;
end

% regional recalls
TP_ISPB       = nnz(pos_IS &  yhat);
TN_OSPB       = nnz(neg_OS & ~yhat);

REC_ISPB_wet  = TP_ISPB / max(1, N_wet_ISPB);
REC_OSPB_dry  = TN_OSPB / max(1, N_dry_OSPB);

REC_ISPB_wet  = max(0,min(1,REC_ISPB_wet));
REC_OSPB_dry  = max(0,min(1,REC_OSPB_dry));

gmean = sqrt(max(REC_ISPB_wet,0) * max(REC_OSPB_dry,0));

if verbose
    fprintf(['[eval:Gmin] n=%d | pos=%d neg=%d | ' ...
             'REC_ISPB_wet=%.3f (N=%d) | REC_OSPB_dry=%.3f (N=%d) | ' ...
             'G_ISPBwet×OSPBdry=%.3f | [global REC=%.3f SPEC=%.3f G=%.3f]\n'], ...
        n, npos, nneg, REC_ISPB_wet, N_wet_ISPB, REC_OSPB_dry, N_dry_OSPB, ...
        gmean, rec_all, spec_all, G_global);
end

% Optional heavy bootstrap (skipped by default)
if ~cfg.skip_bootstrap
    try
        cfg_eval = struct('synthetic_G',double(G), 'synthetic_name','MCMC', ...
            'region',struct('mode','ALL'), 'outdir',cfg.outdir, 'run_id','', ...
            'overwrite',true, 'S_static',S_static, 'skip_static_build',true, ...
            'save_artifacts',false, 'uncertainty',struct('mode','analytic','n_boot',0), ...
            'decision_margin_mode','fixed','decision_margin',0, ...
            'residual',struct('enable',false), ...
            'spec_thresh',0.2,'v_keep',10, ...
            'eval_mask_override',S_static.roi_mask, 'labels_override',S_static.y_raw_vec);
        bootstrap_ghf_component(cfg_eval);
    catch ME
        fprintf('[bootstrap skipped] %s\n', ME.message);
    end
end
end

function [rec, spec] = local_confusion(y, yhat)
TP = nnz( y &  yhat); FN = nnz( y & ~yhat);
TN = nnz(~y & ~yhat); FP = nnz(~y &  yhat);
rec  = TP / max(1, (TP+FN));
spec = TN / max(1, (TN+FP));
rec  = max(0,min(1,rec));
spec = max(0,min(1,spec));
end

function lp = log_prior(params, cfg)
switch lower(cfg.prior_type)
case 'softbox'
    lp = -softbox_penalty(params(1), cfg.mu_range, cfg.soft_k);
    lp = lp - softbox_penalty(params(2), cfg.het_range, cfg.soft_k);
    if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
        lp = lp - softbox_penalty(params(3), cfg.theta_range_deg, cfg.soft_k);
    end
case 'gaussian'
    mu0 = mean(cfg.mu_range);  s_mu = diff(cfg.mu_range)/2;
    h0  = mean(cfg.het_range); s_h  = diff(cfg.het_range)/2;
    lp  = -0.5*((params(1)-mu0)/s_mu).^2 -0.5*((params(2)-h0)/s_h).^2;
    if strcmpi(cfg.template_mode,'gradient') && cfg.sample_angle
        t0 = mean(cfg.theta_range_deg); s_t = diff(cfg.theta_range_deg)/2;
        lp = lp - 0.5*((wrap_centered(params(3),t0)-t0)/s_t).^2;
    end
otherwise
    error('Unknown prior_type: %s', cfg.prior_type);
end

% optional clipping penalty (fast analytic proxy)
if isfield(cfg,'penalize_clip') && cfg.penalize_clip
    frac_clip = 0;
    if strcmpi(cfg.template_mode,'gradient')
        mu = params(1); A = params(2); if cfg.A_nonneg, A = max(A,0); end
        if ~isempty(cfg.floor_G)
            u_floor = (cfg.floor_G - mu) / max(A, eps);
            span = max(0, min(1, u_floor) - (-1));  % portion of U in [-1,1] below floor
            frac_clip = max(frac_clip, min(1, span/2));
        end
        if ~isempty(cfg.cap_G)
            u_cap = (cfg.cap_G - mu) / max(A, eps);
            span_hi = max(0, 1 - max(-1, u_cap));   % portion above cap
            frac_clip = max(frac_clip, min(1, span_hi/2));
        end
    else
        mu = params(1); sig = max(params(2), 1e-6);
        if ~isempty(cfg.floor_G)
            zf = (cfg.floor_G - mu) / sig;  frac_clip = max(frac_clip, 0.5*erfc(zf / sqrt(2)));  % Phi(zf)
        end
        if ~isempty(cfg.cap_G)
            zc = (cfg.cap_G - mu) / sig;    frac_clip = max(frac_clip, 0.5*erfc(-zc / sqrt(2))); % 1-Phi(zc)
        end
    end
    lp = lp - (frac_clip.^2) / max(1e-6, cfg.clip_k);
end
end

function p = softbox_penalty(x, rng, k)
a = rng(1); b = rng(2);
if x<a, p = ((a-x)/k).^2; elseif x>b, p = ((x-b)/k).^2; else, p = 0; end
end

function th = wrap_theta(th, range_deg)
lo = range_deg(1); hi = range_deg(2); span = hi - lo;
th = lo + mod(th - lo, 2*span);
if th>hi, th = hi - (th - hi); end
end

function y = wrap_centered(x, c)
y = x;
while y > c+180, y = y-360; end
while y < c-180, y = y+360; end
end

function z = clamp01(z), z = max(0, min(1, z)); end

function ok = save_one_map(G, X, Y, cfg, outdir, tag, clim, keep_open)
% SAVE_ONE_MAP  Save a map and (optionally) keep the last one open.
% Uses display downsampling only; evaluation math untouched.

    if nargin < 7 || isempty(clim)
        g = G(isfinite(G));
        if isempty(g)
            clim = [0 1];
        else
            p = prctile(g,[2 98]);
            clim = [p(1) p(2)];
            if ~isfinite(clim(1)) || ~isfinite(clim(2)) || clim(1)==clim(2)
                clim = [min(g) max(g)];
                if clim(1)==clim(2), clim = clim + [-1 1]; end
            end
        end
    end
    if nargin < 8, keep_open = false; end

    % --- display-only downsampling (keeps evaluation math intact) ---
    ds = 1;
    if isfield(cfg,'map_downsample') && ~isempty(cfg.map_downsample) && cfg.map_downsample > 1
        ds = round(cfg.map_downsample);
    end
    if ds > 1
        Gd = G(1:ds:end, 1:ds:end);
        Xmin = min(X(:)); Xmax = max(X(:));
        Ymin = min(Y(:)); Ymax = max(Y(:));
        ok = save_map_figure(Gd, [Xmin Xmax], [Ymin Ymax], clim, cfg, outdir, tag, keep_open);
    else
        Xmin = min(X(:)); Xmax = max(X(:));
        Ymin = min(Y(:)); Ymax = max(Y(:));
        ok = save_map_figure(G, [Xmin Xmax], [Ymin Ymax], clim, cfg, outdir, tag, keep_open);
    end
end

function ok = save_map_figure(G, XYXlim, XYYlim, clim, cfg, outdir, name, keep_open)
% SAVE_MAP_FIGURE  Render and save map; optionally keep figure open.
% XYXlim = [xmin xmax], XYYlim = [ymin ymax]. Avoids pcolor/webgl overload.

    if nargin < 8, keep_open = false; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    f = figure('Visible', ternary(keep_open,'on','off'), ...
               'Color','w','Units','pixels','Position',[100 100 1000 800], ...
               'Renderer','painters', 'GraphicsSmoothing','off');

    imagesc('XData', XYXlim, 'YData', XYYlim, 'CData', G);
    set(gca,'YDir','normal'); axis image tight
    caxis(clim);
    cb = colorbar; cb.Label.String = 'Geothermal heat flow (mW m^{-2})';
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('GHF map - %s', string(name)));

    if isfield(cfg,'fontname'), set(gca,'FontName',cfg.fontname); end
    if isfield(cfg,'fontsize'), set(gca,'FontSize',cfg.fontsize); end

    png_path = fullfile(outdir, sprintf('%s.png', string(name)));
    exportgraphics(f, png_path, 'Resolution', 300);
    if isfield(cfg,'export_pdf') && cfg.export_pdf
        pdf_path = fullfile(outdir, sprintf('%s.pdf', string(name)));
        exportgraphics(f, pdf_path, 'ContentType','vector');
    end

    if ~keep_open
        close(f);
    else
        disp(['[save_map_figure] Keeping figure open: ', char(name)]);
        figure(f);
    end
    ok = true;
end

function U = unit_gradient_angle(X, Y, theta_deg)
t  = deg2rad(theta_deg);
nx = cos(t); ny = sin(t);
xc = X - mean(X(:),'omitnan');
yc = Y - mean(Y(:),'omitnan');
P  = xc*nx + yc*ny;
pmin = min(P(:)); pmax = max(P(:));
if ~isfinite(pmin) || ~isfinite(pmax) || pmax==pmin
    U = zeros(size(P),'like',X);
else
    U = single( 2*(P - (pmin+pmax)/2) / (pmax - pmin) ); % range ~ [-1,1]
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = fmtv(v)
s = sprintf('[%s]', strjoin(compose('%.3f', v(:).'), ', '));
end

function Zs = standardize_Z(Z)
m = mean(Z(:), 'omitnan');
s = std(Z(:), 0, 'omitnan');
if ~isfinite(s) || s < eps('single'), Zs = Z - m; else, Zs = (Z - m) ./ s; end
end

function names = param_names(cfg)
if strcmpi(cfg.template_mode,'gradient')
    names = ternary(cfg.sample_angle, {'mu','A','theta_deg'}, {'mu','A'});
else
    names = {'mu','sigma'};
end
end

function plot_trace(y, labelStr, ~, run_id, outdir)
fig = figure('Color','w','Name',['Trace ' labelStr],'Visible','off');
plot(y,'-'); grid on; xlabel('iter'); ylabel(labelStr);
title(['Trace: ' labelStr]);
outpng = fullfile(outdir, sprintf('mcmc_trace_%s_%s.png', regexprep(lower(labelStr),'[^a-z0-9]+','_'), run_id));
exportgraphics(fig, outpng, 'Resolution', 280);
close(fig);
fprintf('[plot] Saved: %s\n', outpng);
end

function plot_posterior(samples, names, ~, run_id, outdir)
fig = figure('Color','w','Name','Posterior','Visible','off');
n = size(samples,2);
tiledlayout(n,1,'Padding','compact','TileSpacing','compact');
for i=1:n
    nexttile; histogram(samples(:,i), 'Normalization','pdf','NumBins', 40);
    grid on; xlabel(names{i}); ylabel('pdf'); title(['Posterior of ' names{i}]);
end
outpng = fullfile(outdir, sprintf('mcmc_posterior_%s.png', run_id));
exportgraphics(fig, outpng, 'Resolution', 280); close(fig);
fprintf('[plot] Saved: %s\n', outpng);
end

function plot_pairs(samples, names, ~, run_id, outdir)
fig = figure('Color','w','Name','Pairs','Visible','off');
n = size(samples,2);
tiledlayout(n, n, 'Padding','compact','TileSpacing','compact');
for i=1:n
    for j=1:n
        nexttile;
        if i==j
            histogram(samples(:,i), 'Normalization','pdf','NumBins', 30);
            ylabel('pdf'); xlabel(names{i});
        else
            scatter(samples(:,j), samples(:,i), 8, 'k', 'filled', ...
                'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.1);
            xlabel(names{j}); ylabel(names{i});
        end
        grid on;
    end
end
outpng = fullfile(outdir, sprintf('mcmc_pairs_%s.png', run_id));
exportgraphics(fig, outpng, 'Resolution', 280); close(fig);
fprintf('[plot] Saved: %s\n', outpng);
end

function [f_floor, f_cap] = physics_clip_stats(G, flo, cap, roi_idx)
if nargin < 4 || isempty(roi_idx)
    mask = isfinite(G);
else
    mask = false(size(G)); mask(roi_idx) = true; mask = mask & isfinite(G);
end
gg = G(mask); n = numel(gg);
if n == 0, f_floor = NaN; f_cap = NaN; return; end
f_floor = nnz(gg <= flo + eps(class(G))) / n;
f_cap   = nnz(gg >= cap - eps(class(G))) / n;
end

