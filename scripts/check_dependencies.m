function check_dependencies(opts)

%CHECK_DEPENDENCIES  Verify required toolboxes and Antarctic Mapping Tools.
%   check_dependencies()                          % defaults
%   check_dependencies(struct('autoInstall',true,'needAMT',true,'verbose',true))
%
%   Toolboxes checked (must-have):
%     - Image Processing Toolbox  (bwconncomp, bwboundaries)
%     - Statistics & Machine Learning Toolbox (glmfit, glmval, zscore, quantile/prctile)
%   Optional (warn if used by your cfg):
%     - Parallel Computing Toolbox (parfor when cfg.use_parallel=true)
%   Add-on (optional but requested): Antarctic Mapping Tools (AMT)

    % ---------- options ----------
    if nargin<1 || ~isstruct(opts), opts = struct; end
    def.autoInstall = true;    % try to auto-install AMT from GitHub zip
    def.needAMT     = true;    % require AMT (set false to only warn)
    def.verbose     = true;
    opts = fill_defaults(opts, def);

    log = @(varargin) fprintf('[deps] %s\n', sprintf(varargin{:}));
    if ~opts.verbose, log = @(varargin) []; end

    % ---------- MATLAB release sanity ----------
    try
        if verLessThan('matlab','9.8')  % R2020a introduced exportgraphics
            log('MATLAB < R2020a detected; exportgraphics may be unavailable. Using print() fallback.');
        end
    catch
    end

    % ---------- must-have toolboxes ----------
    needs = [ ...
        toolbox_req('Image Processing Toolbox',   'image_toolbox',      {@() existsAll({'bwconncomp','bwboundaries'})}); ...
        toolbox_req('Statistics and Machine Learning Toolbox', 'statistics_toolbox', {@() existsAll({'glmfit','glmval','zscore','prctile'})}); ...
    ];

    % Optional toolbox (warn-only)
    optional = [ ...
        toolbox_req('Parallel Computing Toolbox', 'distrib_computing_toolbox', {@() existsAll({'parfor'})}, false) ...
    ];

    % Check required
    for k = 1:numel(needs)
        t = needs(k);
        have = license('test', t.lic) && all(cellfun(@(f) f(), t.tests));
        if ~have
            missing_msg( true, t.name, t.lic, t.tests );
        else
            log('OK: %s', t.name);
        end
    end

    % Check optional
    for k = 1:numel(optional)
        t = optional(k);
        have = license('test', t.lic) && all(cellfun(@(f) f(), t.tests));
        if ~have
            log('Optional: %s not detected. (Only needed if cfg.use_parallel=true)', t.name);
        else
            log('OK (optional): %s', t.name);
        end
    end

    % ---------- Antarctic Mapping Tools (AMT) ----------
    wantAMT = logical(opts.needAMT);
    if wantAMT
        haveAMT = exist('antbounds','file')==2 || exist('bedmap2','file')==2 || exist('measuresps','file')==2;
        if haveAMT
            log('OK: Antarctic Mapping Tools detected.');
        else
            log('AMT not found on path.');
            if opts.autoInstall
                try
                    log('Attempting AMT auto-install (GitHub ZIP)…');
                    tmpzip = fullfile(tempdir, 'Antarctic-Mapping-Tools.zip');
                    url    = 'https://github.com/chadagreene/Antarctic-Mapping-Tools/archive/refs/heads/master.zip';
                    websave(tmpzip, url);
                    target = fullfile(userpath, 'toolboxes', 'AMT');
                    if ~exist(target,'dir'), mkdir(target); end
                    unzip(tmpzip, target);
                    % The repo unzips into a subfolder—add everything:
                    addpath(genpath(target));
                    savepath;
                    haveAMT = exist('antbounds','file')==2 || exist('bedmap2','file')==2;
                    if haveAMT
                        log('AMT installed & added to path.');
                    else
                        error('AMT functions still not found after install.');
                    end
                catch ME
                    if opts.verbose
                        fprintf('AMT auto-install failed');
                    end
                    haveAMT = false;
                end
            end
            if ~haveAMT
                missing_AMT_msg();
            end
        end
    else
        log('AMT check skipped (needAMT=false).');
    end

    % ---------- nested helpers ----------
    function t = toolbox_req(name, lic, tests, must)
        if nargin<4, must = true; end
        t = struct('name',name,'lic',lic,'tests',{tests},'must',must);
    end

    function ok = existsAll(funNames)
        ok = true;
        for ii=1:numel(funNames)
            ok = ok && (exist(funNames{ii},'file')==2);
        end
    end

    function missing_msg(fatal, name, lic, tests)
        feats = {};
        for jj=1:numel(tests)
            % crude feature names from function handles (for user message)
            try
                fstr = func2str(tests{jj});
                feats{end+1} = fstr; %#ok<AGROW>
            catch
            end
        end
        base = sprintf(['Missing required toolbox: %s.\n' ...
                        '  • License id: %s\n' ...
                        '  • Functions checked: bwconncomp/bwboundaries or glmfit/glmval/zscore/prctile\n' ...
                        'Install via: Home > Add-Ons > Get Add-Ons…\n'], name, lic); %#ok<NASGU>
        if fatal
            error('Dependencies:%sMissing', regexprep(name,'\W+',''), ...
                'Missing required toolbox: %s. Please install it via Add-On Explorer.', name);
        else
            warning('Dependencies:%sMissing', regexprep(name,'\W+',''), ...
                'Optional toolbox missing: %s.', name);
        end
    end

    function missing_AMT_msg()
        if wantAMT
            error(['Antarctic Mapping Tools (AMT) not available.\n' ...
                   'Install via: Home > Add-Ons > Get Add-Ons… (search "Antarctic Mapping Tools")\n' ...
                   'OR clone from GitHub and add to path:\n' ...
                   '  https://github.com/chadagreene/Antarctic-Mapping-Tools\n']);
        else
            warning('AMT not found. Some Antarctic plotting utilities may be unavailable.');
        end
    end

    function s = fill_defaults(s, d)
        f = fieldnames(d);
        for i=1:numel(f)
            if ~isfield(s,f{i}) || isempty(s.(f{i}))
                s.(f{i}) = d.(f{i});
            end
        end
    end
end
