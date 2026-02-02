function plot_spec_sweep_all(OUT_or_dir, cfg)
% plot_spec_sweep_all
% Batch-plots your "metric vs model median with CI" figure for EVERY spec_thresh run
% produced by the refactored sweep (bootstrap_ghf_sweep), using your existing
% plot_metrics_vs_modelcval_withCI(S, cfg) function.
%
% INPUTS
%   OUT_or_dir : either
%       (A) OUT struct returned by bootstrap_ghf_sweep()  (preferred), OR
%       (B) a directory that contains spec_* subfolders with artifacts/S_plot_small_*.mat
%
%   cfg : optional plotting cfg override. Common fields:
%       cfg.metrics = {'G_raw','G_ISPBwet_OSPBdry_raw', ...}  % default list below
%       cfg.draw_inner_band = false
%       cfg.label_nonflats  = true/false
%       cfg.exclude_names   = {'MeanNonFlat'}
%       cfg.xlim            = [20 80]
%       cfg.region.mode     = 'ALL'|'ISPB'|'OSPB'|'SPB'
%
% OUTPUTS
%   Saves PNGs into each run's outdir (same behavior as your plot function).
%
% EXAMPLES
%   OUT = bootstrap_ghf_sweep();
%   plot_spec_sweep_all(OUT);
%
%   plot_spec_sweep_all('figs_out');  % if you only have folders + artifacts
%
% NOTE
%   This assumes plot_metrics_vs_modelcval_withCI(S_plot, cfg_plot) is on path.

if nargin < 2 || isempty(cfg), cfg = struct(); end

% defaults
cfg = setdef(cfg, 'metrics', { ...
    'G_ISPBwet_OSPBdry_raw', ...
    'REC_ISPB_wet_raw', ...
    'REC_OSPB_dry_raw', ...
    'G_raw', ...
    'REC_raw', 'SPEC_raw', 'PREC_raw', 'F1_raw', 'ACC_raw' ...
});
cfg = setdef(cfg, 'xlim', [20 80]);
cfg = setdef(cfg, 'verbose', true);

% Resolve runs: a struct array with fields {spec_thresh, outdir, S_plot, cfg_plot}
RUNS = resolve_runs(OUT_or_dir, cfg);

if isempty(RUNS)
    error('plot_spec_sweep_all:NoRunsFound', ...
        'Could not find any runs to plot (no OUT.runs or no artifacts/S_plot_small_*.mat).');
end

% Loop runs
for r = 1:numel(RUNS)
    R = RUNS(r);

    if cfg.verbose
        fprintf('\n[plot-sweep] Run %d/%d | spec_thresh=%.2f | outdir=%s\n', ...
            r, numel(RUNS), R.spec_thresh, R.outdir);
    end

    % Each run has its own cfg_plot/outdir/run_id; we override metric per plot
    base_cfg = R.cfg_plot;
    base_cfg.outdir = R.outdir;           % ensure consistent
    base_cfg.run_id = char(datetime('now','Format','yyyyMMdd_HHmmss'));

    % carry user overrides that your plot function understands
    base_cfg = merge_structs(base_cfg, cfg);

    % enforce xlim inside the plot function? (your plot function hardcodes xlim([20 80]))
    % We'll keep your hardcode, but also store cfg.xlim in case you remove hardcode later.

    for m = 1:numel(cfg.metrics)
        metric = cfg.metrics{m};

        cfg_k = base_cfg;
        cfg_k.metric = metric;

        % Make run_id unique so plots don't overwrite
        cfg_k.run_id = sprintf('%s_spec%0.2f_%s', base_cfg.run_id, R.spec_thresh, metric);

        try
            plot_metrics_vs_modelcval_withCI(R.S_plot, cfg_k);
        catch ME
            warning('[plot-sweep] FAILED spec=%.2f metric=%s : %s', R.spec_thresh, metric, ME.message);
        end
        close(gcf); % keep batch run from accumulating figures
    end
end

fprintf('\n[plot-sweep] Done. Generated plots for %d runs x %d metrics.\n', numel(RUNS), numel(cfg.metrics));

end

%% =======================================================================
%  Resolve runs from either OUT struct or directory
% =======================================================================
function RUNS = resolve_runs(OUT_or_dir, cfg)
RUNS = struct('spec_thresh',{},'outdir',{},'S_plot',{},'cfg_plot',{});

if isstruct(OUT_or_dir) && isfield(OUT_or_dir,'runs') && isfield(OUT_or_dir,'spec_thresh_list')
    OUT = OUT_or_dir;

    % Preferred: if OUT.runs contains results structs with .cfg and (optionally) artifacts already saved,
    % we still load S_plot_small from disk because that is what plot_metrics... expects.
    for k = 1:numel(OUT.runs)
        st = OUT.spec_thresh_list(k);

        % infer outdir
        outdir = '';
        if iscell(OUT.runs) && ~isempty(OUT.runs{k}) && isstruct(OUT.runs{k}) && isfield(OUT.runs{k},'cfg')
            outdir = OUT.runs{k}.cfg.outdir;
        elseif isfield(OUT,'cfg') && isfield(OUT.cfg,'outdir')
            outdir = fullfile(OUT.cfg.outdir, sprintf('spec_%0.2f', st));
        else
            outdir = fullfile('figs_out', sprintf('spec_%0.2f', st));
        end

        [S_plot, cfg_plot] = load_latest_S_plot(outdir);

        if isempty(S_plot)
            warning('[plot-sweep] Missing S_plot_small in %s (spec=%.2f). Skipping.', outdir, st);
            continue;
        end

        RUNS(end+1) = struct('spec_thresh',st,'outdir',outdir,'S_plot',S_plot,'cfg_plot',cfg_plot); %#ok<AGROW>
    end
    return;
end

% Directory mode
if ischar(OUT_or_dir) || isstring(OUT_or_dir)
    base = char(OUT_or_dir);
    if ~exist(base,'dir')
        error('plot_spec_sweep_all:BadDir','Directory not found: %s', base);
    end

    d = dir(fullfile(base, 'spec_*'));
    d = d([d.isdir]);
    if isempty(d)
        % also allow passing a spec_* folder directly
        if startsWith(string(base), "spec_")
            d = dir(base);
            d = d([d.isdir]);
        else
            warning('No spec_* subfolders found under %s', base);
            return;
        end
    end

    % sort by numeric threshold parsed from folder name
    st = nan(numel(d),1);
    for i=1:numel(d)
        st(i) = parse_spec_from_name(d(i).name);
    end
    [~,ord] = sort(st);
    d = d(ord); st = st(ord);

    for i=1:numel(d)
        outdir = fullfile(d(i).folder, d(i).name);
        [S_plot, cfg_plot] = load_latest_S_plot(outdir);
        if isempty(S_plot), continue; end
        RUNS(end+1) = struct('spec_thresh',st(i),'outdir',outdir,'S_plot',S_plot,'cfg_plot',cfg_plot); %#ok<AGROW>
    end
    return;
end

error('plot_spec_sweep_all:BadInput','OUT_or_dir must be OUT struct or directory path.');

end

function st = parse_spec_from_name(name)
% expects "spec_0.10" etc.
st = NaN;
tok = regexp(name, 'spec_([0-9]*\.?[0-9]+)', 'tokens', 'once');
if ~isempty(tok)
    st = str2double(tok{1});
end
end

function [S_plot, cfg_plot] = load_latest_S_plot(outdir)
S_plot = []; cfg_plot = struct();

artdir = fullfile(outdir, 'artifacts');
if ~exist(artdir,'dir'), return; end

ff = dir(fullfile(artdir, 'S_plot_small_*.mat'));
if isempty(ff), return; end

% newest
[~,ix] = max([ff.datenum]);
f = fullfile(ff(ix).folder, ff(ix).name);

L = load(f);
if isfield(L,'S_plot'), S_plot = L.S_plot; end
if isfield(L,'cfg_plot'), cfg_plot = L.cfg_plot; end

% ensure outdir is the run outdir
cfg_plot.outdir = outdir;

end

%% =======================================================================
%  Tiny helpers (local)
% =======================================================================
function cfg = setdef(cfg, field, val)
if ~isfield(cfg,field) || isempty(cfg.(field)), cfg.(field) = val; end
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
