function S = add_g_uncertainty(~, varargin)
% ADD_G_UNCERTAINTY  Compute σ(Gmin) via FD propagation and aggregate σ(GΔ) per model.
% Usage:
%   S = add_g_uncertainty('datasets_for_gmin/gmin_data.mat', ...
%         'sigma_Ts',1.5,'sigma_Mb',0.02,'sigma_H',50, ...
%         'delta_Ts',0.5,'delta_Mb',0.01,'delta_H',5, ...
%         'cache_dir','', 'use_cache',true, 'save',true,'outfile','')

datafile = 'datasets_for_gmin/gmin_data.mat';
addpath('/Users/megankerr/Documents/Location-oldest-ice/ghf/utils');

p = inputParser;
% uncertainty inputs
addParameter(p,'sigma_Ts',1.5,@isnumeric);
addParameter(p,'sigma_Mb',0.02,@isnumeric);
addParameter(p,'sigma_H',50,@isnumeric);
addParameter(p,'delta_Ts',0.5,@isnumeric);
addParameter(p,'delta_Mb',0.01,@isnumeric);
addParameter(p,'delta_H',5,@isnumeric);
% caching & behavior
addParameter(p,'cache_dir','',@(s)ischar(s)||isstring(s));
addParameter(p,'use_cache',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'save',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'outfile','',@(s)ischar(s)||isstring(s));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
opt = p.Results;

logf = @(varargin) ( opt.verbose && fprintf('[add_g_uncertainty] %s\n', sprintf(varargin{:})) );

assert(isfile(datafile), 'Data file not found: %s', datafile);
L = load(datafile,'S'); S = L.S; clear L;

assert(isequal(size(S.Ts), size(S.H)), 'S.Ts size mismatch vs S.H');
assert(isequal(size(S.Mb), size(S.H)), 'S.Mb size mismatch vs S.H');

% ---------- ROI (for consistency with downstream stats) ----------
if isfield(S,'eval_mask') && isequal(size(S.eval_mask), size(S.H))
    Mv = logical(S.eval_mask);
else
    Mv = isfinite(S.H);
end

% ---------- build cfg expected by estimate_sigma_gmin_fd ----------
flat_unc = struct('sigma_Ts',opt.sigma_Ts, 'sigma_Mb',opt.sigma_Mb, 'sigma_H',opt.sigma_H, ...
                  'delta_Ts',opt.delta_Ts, 'delta_Mb',opt.delta_Mb, 'delta_H',opt.delta_H);
cfg_gmin = struct('gmin_unc', flat_unc, ...
                  'eval_mask', Mv, ...
                  'finite_only', true, ...
                  'use_single', isa(S.H,'single'), ...
                  'use_cache', logical(opt.use_cache));

% cache_dir default lives next to the datafile
if isempty(opt.cache_dir)
    data_dir = fileparts(char(datafile)); if isempty(data_dir) || ~isfolder(data_dir), data_dir = pwd; end
    opt.cache_dir = fullfile(data_dir, '_sigma_cache');
end
if ~exist(opt.cache_dir,'dir'), mkdir(opt.cache_dir); end
cfg_gmin.cache_dir = char(opt.cache_dir);
cfg_gmin.cache_key = sprintf('H%dx%d_single%d', size(S.H,1), size(S.H,2), isa(S.H,'single'));

% ---------- compute σ(Gmin) ----------
if exist('estimate_sigma_gmin_fd','file')==2
    [S.sigma_gmin, S.sigma_gmin_scalar] = estimate_sigma_gmin_fd(S, cfg_gmin);
else
    warning('estimate_sigma_gmin_fd.m not on path. Skipping σ(Gmin).');
    S.sigma_gmin = [];
    S.sigma_gmin_scalar = NaN;
end

% ---------- aggregate σ(GΔ) per model ----------
nM = numel(S.names);
validNames = cellfun(@(c) matlab.lang.makeValidName(c), S.names, 'uni', 0);

% optional fallback scalar per-model uncertainties (mW/m^2)
fallback_err = struct();  % e.g., fallback_err.Hazzard = 15;

nGood = nnz(Mv);
S.sigma_gdif_mean   = nan(nM,1);
S.sigma_gdif_median = nan(nM,1);

if ~isempty(S.sigma_gmin) && isequal(size(S.sigma_gmin), size(S.H))
    sigG_vec = double(S.sigma_gmin(Mv));
else
    sigG_vec = zeros(nGood,1);
end

for i = 1:nM
    field_i = validNames{i};
    % per-pixel σ for model, if available
    if isfield(S,'unc') && isfield(S.unc, field_i) && ~isempty(S.unc.(field_i)) ...
            && isequal(size(S.unc.(field_i)), size(S.H))
        tmp = S.unc.(field_i);
        sigM_vec = double(tmp(Mv));
    % derive σ from bounds for Losing
    elseif strcmp(field_i,'Losing') && isfield(S,'bounds') ...
            && all(isfield(S.bounds, {'Losing_min','Losing_max'})) ...
            && isequal(size(S.bounds.Losing_min), size(S.H)) ...
            && isequal(size(S.bounds.Losing_max), size(S.H))
        half95 = 0.5*abs(S.bounds.Losing_max - S.bounds.Losing_min);
        sigM_vec = double(half95(Mv))/1.96;
    else
        scalar = 0;
        if isfield(fallback_err, field_i), scalar = double(fallback_err.(field_i)); end
        sigM_vec = scalar * ones(nGood,1);
    end

    sgd = sqrt(sigM_vec.^2 + sigG_vec.^2);
    S.sigma_gdif_mean(i)   = mean(sgd,'omitnan');
    S.sigma_gdif_median(i) = median(sgd,'omitnan');

    if opt.verbose
        fprintf('[add_g_uncertainty] σ(GΔ) mean %-16s : %.2f mW m^{-2}\n', field_i, S.sigma_gdif_mean(i));
    end
end

% ---------- save ----------
fprintf('Saving uncertainty params');
outfile = char(opt.outfile); if isempty(outfile), outfile = datafile; end
save(outfile, 'S', '-v7.3');
logf('Saved -> %s', outfile);

end
