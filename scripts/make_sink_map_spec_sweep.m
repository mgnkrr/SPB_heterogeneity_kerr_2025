function make_sink_map_spec_sweep(root_outdir, varargin)
% MAKE_SINK_MAP_SPEC_SWEEP  Run make_sink_map() for each spec_* folder.
%
% Example:
%   make_sink_map_spec_sweep('figs_out', ...
%       'map_outdir', 'maps_out', ...
%       'ispb_shp', 'datasets_for_gmin/inner_spb.shp', ...
%       'ospb_shp', 'datasets_for_gmin/outer_spb.shp', ...
%       'residual', true, ...
%       'visible', false, ...
%       'export_pdf', true, ...
%       'pdf_vector', true);
%
% Assumes each folder looks like:
%   figs_out/spec_0.10/artifacts/S_plot_small_YYYYMMDD_HHMMSS.mat

p = inputParser;
addParameter(p,'map_outdir','maps_out',@(s)ischar(s)||isstring(s));
addParameter(p,'pattern','spec_*',@(s)ischar(s)||isstring(s));
addParameter(p,'artifact_glob','S_plot_small_*.mat',@(s)ischar(s)||isstring(s));

% pass-through options to make_sink_map
addParameter(p,'ispb_shp','',@(s)ischar(s)||isstring(s));
addParameter(p,'ospb_shp','',@(s)ischar(s)||isstring(s));
addParameter(p,'residual',true,@islogical);
addParameter(p,'visible',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'fontsize',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'fontname','Times New Roman',@(s)ischar(s)||isstring(s));
addParameter(p,'contours',[1000 1500 2000 2500 3000 3500 4000],@(v)isnumeric(v));
addParameter(p,'dpi',300,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'export_pdf',false,@islogical);
addParameter(p,'pdf_vector',true,@islogical);

parse(p,varargin{:});
opt = p.Results;

root_outdir = char(root_outdir);
assert(exist(root_outdir,'dir')==7, 'Root outdir not found: %s', root_outdir);

% Find spec folders
D = dir(fullfile(root_outdir, char(opt.pattern)));
D = D([D.isdir]);
D = D(~ismember({D.name},{'.','..'}));

if isempty(D)
    fprintf('[sink-sweep] No spec folders found under %s (pattern=%s)\n', root_outdir, opt.pattern);
    return;
end

fprintf('[sink-sweep] Found %d spec folders in %s\n', numel(D), root_outdir);

for k = 1:numel(D)
    spec_name = D(k).name;  % e.g. spec_0.10
    spec_dir  = fullfile(root_outdir, spec_name);
    art_dir   = fullfile(spec_dir, 'artifacts');

    fprintf('\n[sink-sweep] %d/%d | %s\n', k, numel(D), spec_name);

    % newest artifact in that folder
    A = dir(fullfile(art_dir, char(opt.artifact_glob)));
    if isempty(A)
        warning('[sink-sweep] No artifacts found in %s', art_dir);
        continue;
    end
    [~,ix] = max([A.datenum]);
    f_art = fullfile(art_dir, A(ix).name);

    fprintf('[sink-sweep] Loading %s\n', f_art);
    L = load(f_art, 'S_plot', 'cfg_plot');

    if ~isfield(L,'S_plot')
        warning('[sink-sweep] Missing S_plot in %s', f_art);
        continue;
    end

    Splot = L.S_plot;

    % Make per-spec output dir (avoids overwriting & keeps stuff tidy)
    outdir_k = fullfile(char(opt.map_outdir), spec_name);
    if ~exist(outdir_k,'dir'), mkdir(outdir_k); end

    % run_id tag for filenames/shapefile
    run_id = spec_name;
    if isfield(L,'cfg_plot') && isfield(L.cfg_plot,'run_id') && ~isempty(L.cfg_plot.run_id)
        % keep it short but unique if you want:
        % run_id = sprintf('%s_%s', spec_name, L.cfg_plot.run_id);
        run_id = spec_name;
    end

    try
        make_sink_map(Splot, ...
            'outdir', outdir_k, ...
            'run_id', run_id, ...
            'ispb_shp', opt.ispb_shp, ...
            'ospb_shp', opt.ospb_shp, ...
            'residual', opt.residual, ...
            'visible', opt.visible, ...
            'fontsize', opt.fontsize, ...
            'fontname', opt.fontname, ...
            'contours', opt.contours, ...
            'dpi', opt.dpi, ...
            'export_pdf', opt.export_pdf, ...
            'pdf_vector', opt.pdf_vector);

        fprintf('[sink-sweep] OK: %s\n', spec_name);

    catch ME
        warning('[sink-sweep] FAILED %s : %s', spec_name, ME.message);
    end
end

fprintf('\n[sink-sweep] Done.\n');
end
