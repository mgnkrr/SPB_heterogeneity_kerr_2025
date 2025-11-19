function S = summary_best_region(keep, best, varargin)
% SUMMARY_BEST_REGION  CIs around the best G-mean region.
% Usage:
%   L = load('figs_out/artifacts/mcmc_results_YYYYMMDD_HHMMSS.mat');
%   S = summary_best_region(L.keep, L.best, 'mode','topfrac', 'top_frac',0.05);
%   % or nearest-neighbor region:
%   S = summary_best_region(L.keep, L.best, 'mode','knn', 'k_neighbors',500);
%
% Inputs:
%   keep : table from your MCMC (R.keep)
%   best : struct with fields .gmean and .params (R.best)
%
% Options:
%   'mode'         : 'topfrac' (default) or 'knn'
%   'top_frac'     : fraction for topfrac (default 0.05)
%   'k_neighbors'  : integer for knn (default 500)
%   'hpd'          : include HPD95 (default true)
%   'verbose'      : print summary (default true)
%
% Output:
%   S.region.mask    : logical index over keep
%   S.Gmean.(...)    : median, CI50, CI95, HPD95, MCSE, ESS
%   S.params.(name)  : same summaries per parameter in region
%   S.meta           : details (mode, thresholds, sizes, etc.)

% ---------- parse options ----------
p = inputParser;
addParameter(p,'mode','topfrac',@(s)ischar(s)||isstring(s));
addParameter(p,'top_frac',0.05,@(x)isnumeric(x)&&x>0&&x<1);
addParameter(p,'k_neighbors',500,@(x)isnumeric(x)&&x>=10);
addParameter(p,'hpd',true,@islogical);
addParameter(p,'verbose',true,@islogical);
parse(p,varargin{:});
opt = p.Results;

% ---------- checks ----------
assert(istable(keep) && ismember('Gmean', keep.Properties.VariableNames), ...
    'keep must be a table with column Gmean.');
assert(isfield(best,'gmean') && isfield(best,'params'), ...
    'best must contain fields gmean and params.');

% identify parameter columns present
param_cols = intersect({'mu','A','theta_deg','sigma'}, keep.Properties.VariableNames, 'stable');
assert(~isempty(param_cols), 'No parameter columns (mu/A/theta_deg/sigma) found in keep.');

% ---------- choose region ----------
G = keep.Gmean;
N = numel(G);

switch lower(opt.mode)
    case 'topfrac'
        thr = quantile(G, 1 - opt.top_frac);
        mask = (G >= thr);
        % always include the single best draw (closest params to best.params)
        [~, ibest] = min(param_distance(keep, best.params, param_cols));
        mask(ibest) = true;

    case 'knn'
        d = param_distance(keep, best.params, param_cols);
        [~, order] = sort(d, 'ascend');
        k = min(opt.k_neighbors, N);
        mask = false(N,1); mask(order(1:k)) = true;

    otherwise
        error('Unknown mode: %s', opt.mode);
end

% collect region draws
Rgn = keep(mask, :);
n_rgn = height(Rgn);

% ---------- helpers ----------
q = @(x,p) quantile(x(:), p);

    function I = hpd_interval(x, alpha)
        x = sort(x(:));
        n = numel(x); m = floor((1-alpha)*n);
        if m < 1, I = [NaN NaN]; return; end
        spans = x((1+m):n) - x(1:(n-m));
        [~,k] = min(spans);
        I = [x(k) x(k+m)];
    end

    function [mcse, ess] = mcse_batchmeans(x, nbatch)
        x = x(:);
        N = numel(x);
        if nargin<2, nbatch = max(10, floor(sqrt(N))); end
        b = floor(N/nbatch);
        if b < 2, mcse = NaN; ess = NaN; return; end
        x = x(1:(b*nbatch));
        X = reshape(x, b, nbatch);
        bm = mean(X,1);
        mcse = std(bm,1)/sqrt(nbatch);
        ess = var(x,1)/(mcse^2);
    end

    function S1 = summarize_vec(x)
        x = x(:); x = x(isfinite(x));
        S1 = struct('n',numel(x), ...
            'median', median(x), ...
            'CI50',  [q(x,0.25) q(x,0.75)], ...
            'CI95',  [q(x,0.025) q(x,0.975)], ...
            'HPD95', [NaN NaN], ...
            'MCSE',  NaN, 'ESS', NaN);
        if opt.hpd, S1.HPD95 = hpd_interval(x, 0.05); end
        [S1.MCSE, S1.ESS] = mcse_batchmeans(x);
    end

% ---------- summarize G-mean in region ----------
S = struct();
S.region = struct('mode',opt.mode, 'mask',mask, 'n_total',N, 'n_region',n_rgn);

if strcmpi(opt.mode,'topfrac')
    S.region.top_frac = opt.top_frac;
    S.region.threshold = thr;
else
    S.region.k_neighbors = min(opt.k_neighbors, N);
end

S.Gmean = summarize_vec(Rgn.Gmean);

% ---------- summarize parameters in region ----------
S.params = struct();
for c = string(param_cols)
    S.params.(c) = summarize_vec(Rgn.(c));
end

% ---------- optional: extra metrics you saved ----------
if ismember('G_mu_grid', keep.Properties.VariableNames)
    S.G_mu_grid = summarize_vec(Rgn.G_mu_grid);
end
if ismember('G_median_specvalid', keep.Properties.VariableNames)
    S.G_median_specvalid = summarize_vec(Rgn.G_median_specvalid);
end

% ---------- metadata ----------
S.meta = struct();
S.meta.best = struct('gmean',best.gmean, 'params',best.params(:)');

% ---------- print ----------
if opt.verbose
    fprintf('[best-region] mode=%s | region n=%d of %d\n', opt.mode, n_rgn, N);
    if strcmpi(opt.mode,'topfrac')
        fprintf('  threshold (top %.1f%%) Gmean >= %.3f\n', 100*opt.top_frac, S.region.threshold);
    else
        fprintf('  k-nearest draws: k = %d\n', S.region.k_neighbors);
    end
    fprintf('Best Gmean = %.3f @ params = %s\n', best.gmean, vecstr(best.params));

    % Gmean summary
    g = S.Gmean;
    fprintf('Gmean: median %.3f | 50%% [%.3f, %.3f] | 95%% [%.3f, %.3f] | HPD95 ~ [%.3f, %.3f] | MCSE %.4f | ESS ~ %.0f\n', ...
        g.median, g.CI50(1), g.CI50(2), g.CI95(1), g.CI95(2), g.HPD95(1), g.HPD95(2), g.MCSE, g.ESS);

    % Params
    for c = string(param_cols)
        s = S.params.(c);
        fprintf('%s: median %.2f | 50%% [%.2f, %.2f] | 95%% [%.2f, %.2f] | MCSE %.3f | ESS ~ %.0f\n', ...
            c, s.median, s.CI50(1), s.CI50(2), s.CI95(1), s.CI95(2), s.MCSE, s.ESS);
    end

    % Extra metrics if present
    if isfield(S,'G_mu_grid')
        s = S.G_mu_grid;
        fprintf('G_mu_grid: median %.2f | 95%% [%.2f, %.2f]\n', s.median, s.CI95(1), s.CI95(2));
    end
    if isfield(S,'G_median_specvalid')
        s = S.G_median_specvalid;
        fprintf('G_median_specvalid: median %.2f | 95%% [%.2f, %.2f]\n', s.median, s.CI95(1), s.CI95(2));
    end
end
end

% ---------- helper: param distance ----------
function d = param_distance(keep, pbest, param_cols)
P = zeros(height(keep), numel(param_cols));
for i = 1:numel(param_cols)
    P(:,i) = keep.(param_cols{i});
end
pbest = pbest(:)';
if size(P,2) ~= numel(pbest)
    % If best.params omits theta (or includes it), align by overlapping names
    % Try to rebuild best param vector by column names
    map = containers.Map({'mu','A','theta_deg','sigma'}, nan(1,4));
    % naive fill from order of pbest if sizes match; else assume [mu A theta] ordering
    if numel(pbest)==size(P,2)
        pb = pbest;
    else
        % fallback: truncate or pad with NaN to match
        pb = nan(1,size(P,2));
        nmin = min(numel(pbest), size(P,2));
        pb(1:nmin) = pbest(1:nmin);
    end
    pbest = pb;
end
% standardize columns for isotropic distance
Pz = (P - mean(P,1,'omitnan')) ./ std(P,[],1,'omitnan');
pz = (pbest - mean(P,1,'omitnan')) ./ std(P,[],1,'omitnan');
d = sqrt(sum((Pz - pz).^2, 2));
d(~isfinite(d)) = inf;
end

% ---------- helper: print vector nicely ----------
function s = vecstr(v)
s = sprintf('[%s]', strjoin(compose('%.3f', v(:).'), ', '));
end
