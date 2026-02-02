% =========================
% add_g_fields.m
% =========================
function S = add_g_fields(~, varargin)
% ADD_G_FIELDS  Add Gmin, optional advection correction, and per-model GΔ into S.
% Robust to transposed grids and mismatched sizes; safe(ish) for large fields.

p = inputParser;
addParameter(p,'datafile', fullfile('datasets_for_gmin','gmin_data.mat'), @(s)ischar(s)||isstring(s));
addParameter(p,'outfile',  fullfile('datasets_for_gmin','gmin_data.mat'), @(s)ischar(s)||isstring(s));

addParameter(p,'advection_file', fullfile('datasets_for_gmin','longitudinal_advection_maps.mat'), @(s)ischar(s)||isstring(s));
addParameter(p,'advection_var',  'Delta_q_adv', @(s)ischar(s)||isstring(s));
addParameter(p,'advection_units','mWm2', @(s)any(strcmpi(s,{'mWm2','Wm2'})));
addParameter(p,'apply_advection', true, @(x)islogical(x)||ismember(x,[0 1]));

addParameter(p,'apply_colgan', true, @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'colgan_nc', 'colgan_topo_ghf.nc', @(s)ischar(s)||isstring(s));
addParameter(p,'colgan_cache_mat', fullfile('datasets_for_gmin','Delta_Colgan_onTarget.mat'), @(s)ischar(s)||isstring(s));
addParameter(p,'colgan_x_candidates', "x", @(x)isstring(x)||iscellstr(x));
addParameter(p,'colgan_y_candidates', "y", @(x)isstring(x)||iscellstr(x));
addParameter(p,'colgan_D_candidates', "correction", @(x)isstring(x)||iscellstr(x));

addParameter(p,'keep_gdiff_noadv', false, @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'cast_like_H', true, @(x)islogical(x)||ismember(x,[0 1]));

addParameter(p,'gmin_ref_model', 'Hazzard', @(s)ischar(s)||isstring(s)); % reference HF for Gmin computation
addParameter(p,'save', true, @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'verbose', true, @(x)islogical(x)||ismember(x,[0 1]));
parse(p, varargin{:});
opt = p.Results;

logf = @(varargin) ( opt.verbose && fprintf('[add_g_fields] %s\n', sprintf(varargin{:})) );

% ----------- load S -----------
datafile = char(opt.datafile);
assert(isfile(datafile), 'Data file not found: %s', datafile);
L = load(datafile,'S'); S = L.S; clear L;

% ----------- basic sanity -----------
mustHave = {'Xgrid','Ygrid','H','Ts','Mb','names','models'};
for k = 1:numel(mustHave)
    assert(isfield(S, mustHave{k}), 'S.%s missing in %s', mustHave{k}, datafile);
end
assert(isstruct(S.models), 'S.models must be a struct');

% Hard grid contract: every model must match H (or be a simple transpose)
mf = fieldnames(S.models);
for k = 1:numel(mf)
    A = S.models.(mf{k});
    if ~(isequal(size(A), size(S.H)) || isequal(size(A), fliplr(size(S.H))))
        error('Model %s is %s but H is %s (not a transpose). Regrid this model to H grid first.', ...
              mf{k}, mat2str(size(A)), mat2str(size(S.H)));
    end
end

% ----------- harmonize Ts/Mb to match H (transpose-safe) -----------
[S.Ts, did] = align_to_H(S.Ts, S.H, 'S.Ts'); if did, logf('Transposed S.Ts to match S.H'); end
[S.Mb, did] = align_to_H(S.Mb, S.H, 'S.Mb'); if did, logf('Transposed S.Mb to match S.H'); end

% ----------- cast strategy -----------
likeH = 'double';
if opt.cast_like_H, likeH = class(S.H); end
castlike = @(A) cast(A, likeH);

S.H  = castlike(S.H);
S.Ts = castlike(S.Ts);
S.Mb = castlike(S.Mb);

% ----------- Load/interp Colgan Delta_topo (optional) -----------
Delta_topo = [];
if opt.apply_colgan
    cache_mat = char(opt.colgan_cache_mat);

    if isfile(cache_mat)
        logf('Loading cached Delta_topo: %s', cache_mat);
        C = load(cache_mat);
        if isfield(C,'Delta_topo')
            Delta_topo = C.Delta_topo;
        elseif isfield(C,'Delta')
            Delta_topo = C.Delta;
        else
            warning('Colgan cache mat found but no Delta_topo/Delta inside. Ignoring: %s', cache_mat);
        end
        if ~isempty(Delta_topo)
            [Delta_topo, did] = align_to_H(Delta_topo, S.H, 'Delta_topo(cache)');
            if did, logf('Transposed cached Delta_topo to match S.H'); end
        end
    end

    % If no cache, load from netCDF and interpolate onto S.Xgrid/S.Ygrid
    if isempty(Delta_topo)
        if ~isempty(opt.colgan_nc)
            ncfile = char(opt.colgan_nc);
            logf('Reading Colgan netCDF and interpolating -> S grid: %s', ncfile);

            Delta_topo = read_colgan_and_interp_to_target( ...
                ncfile, S.Xgrid, S.Ygrid, ...
                string(opt.colgan_x_candidates), string(opt.colgan_y_candidates), string(opt.colgan_D_candidates));

            [Delta_topo, did] = align_to_H(Delta_topo, S.H, 'Delta_topo(nc interp)');
            if did, logf('Transposed interpolated Delta_topo to match S.H'); end

            % Save cache for next time (optional but helpful)
            try
                Delta_topo_to_save = Delta_topo; 
                save(cache_mat, 'Delta_topo_to_save', '-v7.3');
                logf('Cached Delta_topo -> %s', cache_mat);
            catch ME
                warning('Could not cache Delta_topo to %s (%s)', cache_mat, ME.message);
            end
        else
            logf('apply_colgan=true but no colgan_nc provided and no cache found. Skipping Colgan.');
        end
    end
end

% store for provenance / later use
S.Delta_topo = castlike(Delta_topo);

if ~isempty(Delta_topo)
    logf('Delta_topo ready: finite=%d/%d | min=%.3g max=%.3g med=%.3g', ...
        nnz(isfinite(Delta_topo)), numel(Delta_topo), ...
        min(Delta_topo(:),[],'omitnan'), max(Delta_topo(:),[],'omitnan'), median(Delta_topo(:),'omitnan'));
end

% ----------- model name normalization & mapping -----------
% Use S.names to define a stable order, but ensure fields exist in S.models.
validNames = cellfun(@matlab.lang.makeValidName, S.names, 'uni', 0);
validNames = matlab.lang.makeUniqueStrings(validNames);

modelFields = fieldnames(S.models);

% Build a mapping from validNames{i} -> actual field in S.models (possibly renamed)
% Strategy:
%   1) exact match on valid name
%   2) case-insensitive match on valid name
%   3) try case-insensitive match on raw original name made valid
% If found as some other field, rename to the valid name for consistency.
for i = 1:numel(validNames)
    vn = validNames{i};
    if isfield(S.models, vn), continue; end

    hit = find(strcmp(modelFields, vn), 1);
    if isempty(hit), hit = find(strcmpi(modelFields, vn), 1); end

    if isempty(hit)
        raw = matlab.lang.makeValidName(S.names{i});
        hit = find(strcmpi(modelFields, raw), 1);
    end

    if isempty(hit)
        error('Model "%s" -> "%s" not found in S.models. Available fields include: %s', ...
            S.names{i}, vn, strjoin(modelFields(1:min(12,end)),', '));
    end

    old = modelFields{hit};
    S.models.(vn) = S.models.(old);
    if ~strcmp(old, vn)
        S.models = rmfield(S.models, old);
        % refresh list after rename
        modelFields = fieldnames(S.models);
    end
end

% ----------- Load & align advection (optional) -----------
D = [];
if opt.apply_advection
    advfile = char(opt.advection_file);
    if ~isempty(advfile) && isfile(advfile)
        logf('Loading advection map: %s (%s)', advfile, opt.advection_var);
        A = load(advfile);
        if isfield(A, opt.advection_var)
            D0 = A.(opt.advection_var);
            if strcmpi(opt.advection_units,'Wm2'), D0 = 1000 * D0; end % -> mW/m^2

            % Try to regrid if source X/Y grids exist
            if all(isfield(A, {'Xgrid','Ygrid'})) && isequal(size(A.Xgrid), size(D0)) && isequal(size(A.Ygrid), size(D0))
                D = regrid_to(S.Xgrid, S.Ygrid, A.Xgrid, A.Ygrid, D0);
            else
                % Otherwise, try direct match / transpose
                [D, ~] = align_to_H(D0, S.H, 'Delta_q_adv');
            end

            D = castlike(D);
            D(~isfinite(D)) = 0; % safe fill (your choice)
        else
            warning('Advection variable "%s" not found in %s. Skipping advection.', opt.advection_var, advfile);
        end
    else
        logf('No advection file found — proceeding without advection.');
    end
else
    logf('apply_advection=false — proceeding without advection.');
end

S.Delta_q_adv = D;

% ----------- Compute Gmin (no-advection) -----------
refName = matlab.lang.makeValidName(char(opt.gmin_ref_model));
refName = matlab.lang.makeUniqueStrings({refName}); refName = refName{1};

% If the ref model name doesn't exist, fall back to first model in validNames.
if ~isfield(S.models, refName)
    warning('Reference model "%s" not found. Using first model "%s" as reference.', refName, validNames{1});
    refName = validNames{1};
end

HFref = S.models.(refName);
[HFref, did] = align_to_H(HFref, S.H, sprintf('S.models.%s', refName));
if did, logf('Transposed reference model "%s" to match S.H', refName); end
HFref = castlike(HFref);

logf('Computing Gmin (no-advection) using reference model: %s', refName);
[~, Gmin_noadv, ~, ~, ~, ~] = processGHF(HFref, S.Ts, S.Mb, S.H, [], S.Delta_topo);

S.Gmin = castlike(Gmin_noadv);

% Also store advected Gmin if D exists
S.Gmin_adv = [];
if ~isempty(D)
    S.Gmin_adv = castlike(S.Gmin - D);
end

% ----------- Per-model GΔ fields -----------
nM = numel(validNames);
if ~isfield(S,'Gdiff') || ~isstruct(S.Gdiff), S.Gdiff = struct(); end
if opt.keep_gdiff_noadv && (~isfield(S,'Gdiff_noadv') || ~isstruct(S.Gdiff_noadv))
    S.Gdiff_noadv = struct();
end

logf('Building GΔ fields for %d models%s ...', nM, tern(~isempty(D),' (with advection)',''));
for i = 1:nM
    field_i = validNames{i};

    Mi = S.models.(field_i);
    [Mi, did] = align_to_H(Mi, S.H, sprintf('S.models.%s', field_i));
    if did && opt.verbose
        fprintf('[add_g_fields] note: transposed model "%s" to match S.H\n', field_i);
    end
    Mi = castlike(Mi);

    if isempty(D)
        [~, ~, Gdiff, ~, ~, ~] = processGHF(Mi, S.Ts, S.Mb, S.H, [], S.Delta_topo);
        S.Gdiff.(field_i) = castlike(Gdiff);
    else
        % Pass D into processGHF (no internal file loads)
        [~, ~, Gdiff, ~, Gdiff_adv, ~] = processGHF(Mi, S.Ts, S.Mb, S.H, D, S.Delta_topo);

        if opt.keep_gdiff_noadv
            S.Gdiff_noadv.(field_i) = castlike(Gdiff);
        end
        S.Gdiff.(field_i) = castlike(Gdiff_adv);
    end

    if opt.verbose
        fprintf('[add_g_fields] %02d/%02d  %-28s ok\n', i, nM, field_i);
    end
end

% ----------- Save (robust) -----------
outfile = char(opt.outfile);
if isempty(outfile), outfile = datafile; end

if opt.save
    logf('Saving to %s ...', outfile);
    try
        mf = matfile(outfile,'Writable',true);
        mf.S = S;  % single assignment; avoids extra copies
        logf('Saved -> %s', outfile);
    catch ME
        fprintf('matfile write failed . Falling back to save -v7.3');
        save(outfile, 'S', '-v7.3');
        logf('Saved (fallback) -> %s', outfile);
    end
end
end

% ---------- helpers ----------
function [Aout, didTranspose] = align_to_H(Ain, H, name)
% Ensure Ain matches size(H). If Ain is transpose, fix it.
didTranspose = false;

if isequal(size(Ain), size(H))
    Aout = Ain;
    return;
end

if isequal(size(Ain), fliplr(size(H)))
    Aout = Ain.';
    didTranspose = true;
    return;
end

error('align_to_H:SizeMismatch', ...
    '%s is %s but H is %s. Provide data on H grid or regrid explicitly.', ...
    name, mat2str(size(Ain)), mat2str(size(H)));
end

function Dout = regrid_to(Xt, Yt, Xs, Ys, Din)
% Robust regridding with monotonic axes
D = double(Din);
Fx = double(Xs(1,:));
Fy = double(Ys(:,1));

% Ensure monotonic increasing axes for griddedInterpolant
if any(diff(Fx) <= 0)
    Fx = fliplr(Fx);
    D  = fliplr(D);
end
if any(diff(Fy) <= 0)
    Fy = flipud(Fy);
    D  = flipud(D);
end

F = griddedInterpolant({Fy, Fx}, D, 'linear', 'nearest'); % {row(Y), col(X)}
Dout = F(double(Yt), double(Xt));
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

function Delta_out = read_colgan_and_interp_to_target(ncfile, Xgrid, Ygrid, x_cand, y_cand, D_cand)
assert(isfile(ncfile), 'Colgan netCDF not found: %s', ncfile);

info = ncinfo(ncfile);
vars = string({info.Variables.Name});

x_name = first_present(vars, x_cand);
y_name = first_present(vars, y_cand);
D_name = first_present(vars, D_cand);

if x_name == "" || y_name == "" || D_name == ""
    fprintf('\n[Colgan] Could not auto-detect required vars in: %s\n', ncfile);
    fprintf('Available variables:\n');
    disp(vars.');
    error('add_g_fields:ColganVarDetect', ...
        'Set colgan_*_candidates to match the netCDF variable names.');
end

x = ncread(ncfile, char(x_name));
y = ncread(ncfile, char(y_name));
D = ncread(ncfile, char(D_name));

x = double(x(:)');   % 1×nc
y = double(y(:));    % nr×1
D = double(D);

% Make D be [numel(y) × numel(x)] (rows=y, cols=x)
if isequal(size(D), [numel(x) numel(y)])
    D = D.'; % transpose to [y×x]
end
assert(isequal(size(D), [numel(y) numel(x)]), ...
    'Colgan Delta size (%dx%d) does not match x/y lengths (%d,%d).', ...
    size(D,1), size(D,2), numel(y), numel(x));

% Ensure x ascending, y ascending
if x(1) > x(end)
    x = fliplr(x);
    D = fliplr(D);
end
if y(1) > y(end)
    y = flipud(y);
    D = flipud(D);
end

% Interpolate onto target
Delta_out = interp2(x, y, D, Xgrid, Ygrid, 'linear', NaN);

fprintf('[Colgan] loaded %s (x=%s y=%s D=%s) -> interpolated to %dx%d\n', ...
    ncfile, x_name, y_name, D_name, size(Xgrid,1), size(Xgrid,2));
end

function name = first_present(allVars, candidates)
name = "";
candidates = string(candidates);
for k = 1:numel(candidates)
    if any(strcmpi(allVars, candidates(k)))
        idx = find(strcmpi(allVars, candidates(k)), 1, 'first');
        name = allVars(idx);
        return;
    end
end
end
