function T = scrape_threshold_table_simple(sweep_out, spec_list)
% Scrape numeric G_raw values (EXP + optional CI) for each model at each threshold.
% - Uses artifacts/S_plot_small_*.mat in the SAME ORDER as spec_list
% - Keeps NON-FLAT models only
%
% Output columns:
%   spec_thresh, model_name, model_title, G_raw, G_raw_lo95, G_raw_hi95, G_raw_lo50, G_raw_hi50
%
% Usage:
%   spec_list = [0.1 0.15 0.2 0.25 0.3 0.35];
%   T = scrape_threshold_table_simple_Graw(sweep_out, spec_list);

if nargin < 1 || isempty(sweep_out), sweep_out = pwd; end
if nargin < 2 || isempty(spec_list), error('Provide spec_list (threshold order).'); end

A = dir(fullfile(sweep_out, 'artifacts', 'S_plot_small_*.mat'));
assert(~isempty(A), 'No S_plot_small_*.mat found under %s/artifacts', sweep_out);

% Sort by filename (timestamp in name) so order is stable
[~,ord] = sort({A.name});
A = A(ord);

assert(numel(A) >= numel(spec_list), ...
    'Found %d artifacts but spec_list has %d thresholds. Need >=.', numel(A), numel(spec_list));

T = table();

for i = 1:numel(spec_list)
    spec  = spec_list(i);
    fpath = fullfile(A(i).folder, A(i).name);

    L = load(fpath, 'S_plot');
    assert(isfield(L,'S_plot') && isstruct(L.S_plot), 'Bad artifact: %s', fpath);
    S = L.S_plot;

    assert(isfield(S,'Gstats') && isfield(S.Gstats,'G_raw'), ...
        'G_raw not found in S_plot.Gstats for %s', fpath);

    [G, lo95, hi95, lo50, hi50, ok] = pull_ci(S.Gstats, 'G_raw');
    assert(ok, 'G_raw missing EXP in %s', fpath);

    names  = string(S.names(:));
    titles = names;
    if isfield(S,'titles') && numel(S.titles)==numel(S.names)
        titles = string(S.titles(:));
    end

    n = min([numel(names), numel(G)]);
    names  = names(1:n);
    titles = titles(1:n);
    G      = G(1:n);
    lo95   = lo95(1:n);
    hi95   = hi95(1:n);
    lo50   = lo50(1:n);
    hi50   = hi50(1:n);

    is_flat = startsWith(names, "Flat_");
    keep    = ~is_flat;

    Ti = table();
    Ti.spec_thresh = repmat(spec, nnz(keep), 1);
    Ti.model_name  = names(keep);
    Ti.model_title = titles(keep);

    % --- NUMERIC G_raw values you care about ---
    Ti.G_raw      = G(keep);
    Ti.G_raw_lo95 = lo95(keep);
    Ti.G_raw_hi95 = hi95(keep);
    Ti.G_raw_lo50 = lo50(keep);
    Ti.G_raw_hi50 = hi50(keep);

    Ti.source_file = repmat(string(A(i).name), nnz(keep), 1);

    T = [T; Ti]; %#ok<AGROW>
end

T = sortrows(T, {'spec_thresh','model_name'});

out_csv = fullfile(sweep_out, sprintf('SCRAPED_SIMPLE_G_raw_%s.csv', char(datetime('now','Format','yyyyMMdd_HHmmss'))));
writetable(T, out_csv);
fprintf('[simple scrape] wrote: %s\n', out_csv);

end

function [Y, lo, hi, lo50, hi50, ok] = pull_ci(Gstats, metric)
Y=[]; lo=[]; hi=[]; lo50=[]; hi50=[]; ok=false;
if ~isfield(Gstats, metric), return; end
Gf = Gstats.(metric);

if ~isfield(Gf,'EXP') || isempty(Gf.EXP), return; end
Y = double(Gf.EXP(:));

lo = NaN(size(Y)); hi = NaN(size(Y));
if isfield(Gf,'CI') && ~isempty(Gf.CI) && size(Gf.CI,2) >= 2
    lo = double(Gf.CI(:,1));
    hi = double(Gf.CI(:,2));
end

lo50 = NaN(size(Y)); hi50 = NaN(size(Y));
if isfield(Gf,'CIi') && ~isempty(Gf.CIi) && size(Gf.CIi,2) >= 2
    lo50 = double(Gf.CIi(:,1));
    hi50 = double(Gf.CIi(:,2));
end

ok = true;
end
