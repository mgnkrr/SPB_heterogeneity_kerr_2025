function S = build_data_core(varargin)
% BUILD_DATA_CORE  Load base datasets into struct S (no Gmin/GΔ).
%
%   S = build_data_core('to_single',true,'save',true,...
%                       'outdir','datasets_for_gmin','outfile','gmin_data.mat')

%% -------- options & logger --------
p = inputParser;
addParameter(p,'to_single',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'save',true,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'outdir','datasets_for_gmin',@(s)ischar(s)||isstring(s));
addParameter(p,'outfile','gmin_data.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'verbose',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p,varargin{:});
opt = p.Results;

logf = @(varargin) ( opt.verbose && fprintf('[build_data_core] %s\n', sprintf(varargin{:})) );

S = struct();

%% -------- GHF models --------
t = tic; logf('Loading GHF models ...');
t1=load('datasets_for_gmin/GHF_Hazzard_interp_corr.mat','GHF_Hazzard_interp_corr');   Hazz    = t1.GHF_Hazzard_interp_corr;    clear t1
t1=load('datasets_for_gmin/GHF_Martos_interp_corr.mat','GHF_Martos_interp_corr');     Mart    = t1.GHF_Martos_interp_corr;     clear t1
t1=load('datasets_for_gmin/GHF_Shen_interp_corr.mat','GHF_Shen_interp_corr');         Shen    = t1.GHF_Shen_interp_corr;       clear t1
t1=load('datasets_for_gmin/GHF_Stal_interp_corr.mat','GHF_Stal_interp_corr');         Stal    = t1.GHF_Stal_interp_corr;       clear t1
t1=load('datasets_for_gmin/GHF_Losing_interp_corr.mat','GHF_Losing_interp_corr');     Los     = t1.GHF_Losing_interp_corr;     clear t1
t1=load('datasets_for_gmin/GHF_An_interp_corr.mat','GHF_An_interp_corr');             An      = t1.GHF_An_interp_corr;         clear t1
t1=load('datasets_for_gmin/GHF_FoxMaule_interp_corr.mat','GHF_FoxMaule_interp_corr'); Fox     = t1.GHF_FoxMaule_interp_corr;   clear t1
t1=load('datasets_for_gmin/GHF_Haeger_interp_corr.mat','GHF_Haeger_interp_corr');     Haeger  = t1.GHF_Haeger_interp_corr;     clear t1
t1=load('datasets_for_gmin/GHF_Lucazeau_interp_corr.mat','GHF_Lucazeau_interp_corr'); Lucazeau = t1.GHF_Lucazeau_interp_corr;  clear t1 

S.models = struct( ...
    'FoxMaule', Fox, ...
    'Hazzard',  Hazz, ...
    'Martos',   Mart, ...
    'Shen',     Shen, ...
    'Stal',     Stal, ...
    'An',       An, ...
    'Losing',   Los, ...
    'Haeger',   Haeger, ...
    'Lucazeau', Lucazeau);

S.names  = {'FoxMaule','Hazzard','Martos','Shen','Stal','An','Losing','Haeger','Lucazeau'};
S.titles = {'Fox Maule','Hazzard','Martos','Shen','Stal','An','Losing','Haeger','Lucazeau'};
logf('GHF models loaded in %.2f s', toc(t));

%% -------- Grid / ice thickness --------
t = tic; logf('Loading grid + H ...');
g = load('datasets_for_gmin/coldex_icethk.mat','Xgrid','Ygrid','H');
S.Xgrid = g.Xgrid; S.Ygrid = g.Ygrid; S.H = g.H; clear g
S.X_km = S.Xgrid/1000; S.Y_km = S.Ygrid/1000;
logf('Grid size: %dx%d | done in %.2f s', size(S.Xgrid,1), size(S.Xgrid,2), toc(t));

%% -------- Specularity (9999 -> NaN) --------
t = tic; logf('Loading specularity ...');
q = load('datasets_for_gmin/specularity.mat','Q'); 
S.Q = q.Q; clear q
S.spec_invalid = ~isfinite(S.Q);
nInvalid = nnz(S.spec_invalid);

logf('Spec invalid (==9999): %d cells (%.2f%%) | %.2f s', nInvalid, 100*nInvalid/max(1,numel(S.Q)), toc(t));

%% -------- Ts, Mb, ice velocity --------
t = tic; logf('Loading Ts, Mb, ice velocity ...');
tmp = load('datasets_for_gmin/Ts_interp.mat','Ts_interp','Xgrid','Ygrid');
S.Ts = tmp.Ts_interp;
assert(isequal(size(S.Ts), size(S.H)), 'Ts size mismatch');
if isfield(tmp,'Xgrid')
    assert(isequal(tmp.Xgrid, S.Xgrid) && isequal(tmp.Ygrid, S.Ygrid), 'Ts Xgrid/Ygrid mismatch');
end
clear tmp

tmp = load('datasets_for_gmin/Mb_interp.mat','Mb_interp','Xgrid','Ygrid');
S.Mb = tmp.Mb_interp;
assert(isequal(size(S.Mb), size(S.H)), 'Mb size mismatch');
if isfield(tmp,'Xgrid')
    assert(isequal(tmp.Xgrid, S.Xgrid) && isequal(tmp.Ygrid, S.Ygrid), 'Mb Xgrid/Ygrid mismatch');
end
clear tmp

tmp = load('datasets_for_gmin/mouginot_icevel.mat','speed','Xgrid','Ygrid');
S.icevel = tmp.speed;
assert(isequal(size(S.icevel), size(S.H)), 'icevel size mismatch');
if isfield(tmp,'Xgrid')
    assert(isequal(tmp.Xgrid, S.Xgrid) && isequal(tmp.Ygrid, S.Ygrid), 'icevel Xgrid/Ygrid mismatch');
end
clear tmp

logf('Ts/Mb/icevel loaded in %.2f s', toc(t));

%% -------- Uncertainty / bounds (optional files) --------
t = tic; logf('Loading uncertainty/bounds if present ...');
S.unc    = struct(); 
S.bounds = struct();

if isfile('datasets_for_gmin/UNC_Hazzard_interp.mat')
    u = load('datasets_for_gmin/UNC_Hazzard_interp.mat','UNC_Hazzard_interp');
    S.unc.Hazzard = u.UNC_Hazzard_interp; clear u;
end
if isfile('datasets_for_gmin/UNC_Martos_interp.mat')
    u = load('datasets_for_gmin/UNC_Martos_interp.mat','UNC_Martos_interp');
    S.unc.Martos = u.UNC_Martos_interp; clear u;
end
if isfile('datasets_for_gmin/UNC_Shen_interp.mat')
    u = load('datasets_for_gmin/UNC_Shen_interp.mat','UNC_Shen_interp');
    S.unc.Shen = u.UNC_Shen_interp; clear u;
end
if isfile('datasets_for_gmin/UNC_Stal_interp.mat')
    u = load('datasets_for_gmin/UNC_Stal_interp.mat','UNC_Stal_interp');
    S.unc.Stal = u.UNC_Stal_interp; clear u;
end
if isfile('datasets_for_gmin/UNC_Haeger_interp.mat')
    u = load('datasets_for_gmin/UNC_Haeger_interp.mat','UNC_Haeger_interp');
    S.unc.Haeger = u.UNC_Haeger_interp; clear u;
end
if isfile('datasets_for_gmin/BMIN_Losing_interp.mat')
    b = load('datasets_for_gmin/BMIN_Losing_interp.mat','BMIN_Losing_interp');
    S.bounds.Losing_min = b.BMIN_Losing_interp; clear b;
end
if isfile('datasets_for_gmin/BMAX_Losing_interp.mat')
    b = load('datasets_for_gmin/BMAX_Losing_interp.mat','BMAX_Losing_interp');
    S.bounds.Losing_max = b.BMAX_Losing_interp; clear b;
end
if isfile('datasets_for_gmin/UNC_Lucazeau_interp.mat')
    u = load('datasets_for_gmin/UNC_Lucazeau_interp.mat','UNC_Lucazeau_interp');
    S.unc.Lucazeau = u.UNC_Lucazeau_interp; clear u;
end

unc_fields = {'Hazzard','Martos','Shen','Stal','Haeger','Lucazeau'};
have_unc = unc_fields(isfield(S.unc, unc_fields));
if ~isempty(have_unc)
    Usum = 0; Ucnt = 0;
    for k = 1:numel(have_unc)
        U = S.unc.(have_unc{k});
        if ~isempty(U)
            Usum = Usum + double(U);
            Ucnt = Ucnt + 1;
        end
    end
end
logf('Uncertainty/bounds load in %.2f s', toc(t));

%% -------- Optional cast to single --------
if opt.to_single
    t = tic; logf('Casting numeric fields to single ...');
    S.H      = single(S.H);
    S.Q      = single(S.Q);
    S.icevel = single(S.icevel);
    S.Ts     = single(S.Ts);
    S.Mb     = single(S.Mb);
    f = fieldnames(S.models);
    for k = 1:numel(f)
        S.models.(f{k}) = single(S.models.(f{k}));
    end
    if isfield(S,'unc')
        f = fieldnames(S.unc);
        for k = 1:numel(f)
            if ~isempty(S.unc.(f{k}))
                S.unc.(f{k}) = single(S.unc.(f{k}));
            end
        end
    end
    if isfield(S,'bounds')
        if isfield(S.bounds,'Losing_min') && ~isempty(S.bounds.Losing_min)
            S.bounds.Losing_min = single(S.bounds.Losing_min);
        end
        if isfield(S.bounds,'Losing_max') && ~isempty(S.bounds.Losing_max)
            S.bounds.Losing_max = single(S.bounds.Losing_max);
        end
    end
    logf('Cast complete in %.2f s', toc(t));
end

%% -------- Flat models --------
logf('Adding flat GHF layers ...');
flat_vals = 20:1:80;  % mW/m^2
for v = flat_vals
    fld = matlab.lang.makeValidName(sprintf('Flat_%0.4g', v));
    F = zeros(size(S.H), 'like', S.H) + cast(v, 'like', S.H);
    S.models.(fld) = F;
    S.names{end+1}  = fld;
    S.titles{end+1} = sprintf('Flat %0.4g', v);
    if ~isfield(S,'unc'), S.unc = struct(); end
    S.unc.(fld) = zeros(size(S.H), 'like', S.H);
end

%% -------- Save (optional) --------
if opt.save
    if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end
    outfile = fullfile(opt.outdir, char(opt.outfile));
    save(outfile, 'S', '-v7.3');
    logf('Saved -> %s', outfile);
end
end
