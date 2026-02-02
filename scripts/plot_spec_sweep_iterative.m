function plot_spec_sweep_iterative(base_outdir, metrics)
% Loop over spec_* folders and call plot_metrics_vs_modelcval_withCI once per spec run.
%
% Assumes each run saved:
%   <base_outdir>/spec_0.10/artifacts/S_plot_small_*.mat
% and that file contains variables: S_plot and cfg_plot
%
% Usage:
%   plot_spec_sweep_iterative('figs_out');
%   plot_spec_sweep_iterative('figs_out', {'G_raw','REC_raw'});

if nargin < 2 || isempty(metrics)
    metrics = {'G_raw'};  % keep it minimal; add more if you want
end

d = dir(fullfile(base_outdir, 'spec_*'));
d = d([d.isdir]);

if isempty(d)
    error('No spec_* folders found in: %s', base_outdir);
end

fprintf('[spec-sweep] Found %d spec folders in %s\n', numel(d), base_outdir);

for ii = 1:numel(d)
    spec_dir = fullfile(base_outdir, d(ii).name);
    art_dir  = fullfile(spec_dir, 'artifacts');

    mats = dir(fullfile(art_dir, 'S_plot_small_*.mat'));
    if isempty(mats)
        warning('[spec-sweep] No S_plot_small_*.mat found in %s', art_dir);
        continue
    end

    % pick newest artifact
    [~, ix] = max([mats.datenum]);
    mat_path = fullfile(art_dir, mats(ix).name);

    fprintf('\n[spec-sweep] %d/%d | %s\n', ii, numel(d), d(ii).name);
    fprintf('[spec-sweep] Loading %s\n', mat_path);

    L = load(mat_path, 'S_plot', 'cfg_plot');

    % Ensure outputs go into the spec_* folder you're iterating over
    cfg = L.cfg_plot;
    cfg.outdir = spec_dir;

    for m = 1:numel(metrics)
        cfg.metric = metrics{m};

        % keep run_id stable per spec folder but unique per metric
        cfg.run_id = sprintf('%s_%s', d(ii).name, cfg.metric);

        try
            plot_metrics_vs_modelcval_withCI(L.S_plot, cfg);
            %close(gcf);
        catch ME
            warning('[spec-sweep] FAILED %s metric=%s : %s', d(ii).name, cfg.metric, ME.message);
        end
    end
end
end
