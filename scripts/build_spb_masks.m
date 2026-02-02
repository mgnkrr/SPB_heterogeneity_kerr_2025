function build_spb_masks()
% BUILD_SPB_MASKS (verbose)
% - Crops to rectangle before any inpolygon
% - Logs timings, sizes, counts
% - Saves either cropped or full masks (cropped by default, much smaller)

%% ------------------ Verbose helpers ------------------
t0 = tic;
memMB = @() struct2table(whos, 'AsArray',true).bytes'/1024/1024;
function logf(fmt, varargin), fprintf(['[build_spb_masks] ' fmt '\n'], varargin{:}); end

logf('Start');

%% ------------------ Config ------------------
cfg.rect = struct('enable',true,'xlim',[-200 800],'ylim',[-200 400],'rects',[],'mode','include');
useCroppedSave = true;    % save only cropped masks + indices (smaller, faster)
sparseThreshold = 0.10;   % save as sparse if nnz < 10%

%% ------------------ Grid ------------------
t = tic;
G = load('datasets_for_gmin/coldex_icethk.mat','Xgrid','Ygrid','H');
Xgrid = G.Xgrid; Ygrid = G.Ygrid; % H unused for masks; keep memory light
clear G
X_km = Xgrid/1000; Y_km = Ygrid/1000;
logf('Loaded grid: %dx%d (%.1f MB)', size(X_km,1), size(X_km,2), sum(memMB()));
logf('Grid load time: %.2f s', toc(t));

%% ------------------ Rectangle mask + crop ------------------
t = tic;
rect_mask = build_rect_mask(X_km, Y_km, cfg.rect);
[rowIdx,colIdx] = find(rect_mask);
if isempty(rowIdx) || isempty(colIdx)
    logf('Rectangle empty; falling back to full domain.');
    rmin=1; rmax=size(X_km,1); cmin=1; cmax=size(X_km,2);
else
    rmin = min(rowIdx); rmax = max(rowIdx);
    cmin = min(colIdx); cmax = max(colIdx);
end
Xc = X_km(rmin:rmax, cmin:cmax);
Yc = Y_km(rmin:rmax, cmin:cmax);
logf('Rect crop rows [%d..%d], cols [%d..%d] -> %dx%d window', rmin,rmax,cmin,cmax,size(Xc,1),size(Xc,2));
logf('Rect build time: %.2f s', toc(t));

%% ------------------ Read shapefiles ------------------
t = tic;
S_spb  = shaperead('datasets_for_gmin/spb_merged.shp');
S_ispb = shaperead('datasets_for_gmin/inner_spb.shp');
S_ospb = shaperead('datasets_for_gmin/outer_spb.shp');
log_parts = @(S) sum(arrayfun(@(k) numel(S(k).X), 1:numel(S)));
logf('Loaded shapes: SPB parts=%d verts=%d | ISPB parts=%d verts=%d | OSPB parts=%d verts=%d', ...
     numel(S_spb),  log_parts(S_spb), ...
     numel(S_ispb), log_parts(S_ispb), ...
     numel(S_ospb), log_parts(S_ospb));
logf('Shapefile load time: %.2f s', toc(t));

%% ------------------ Build masks (on CROPPED grid) ------------------
% We'll OR-accumulate per part. (You can switch to concat->single call if desired.)
mask_SPB_c  = false(size(Xc));
mask_ISPB_c = false(size(Xc));
mask_OSPB_c = false(size(Xc));

% SPB
t = tic;
for i = 1:numel(S_spb)
    in = inpolygon(Xc, Yc, S_spb(i).X/1000, S_spb(i).Y/1000);
    mask_SPB_c = mask_SPB_c | in;
    if mod(i,50)==0 || i==numel(S_spb), logf('SPB %d/%d parts...', i, numel(S_spb)); end
end
logf('SPB mask crop nnz=%d (%.3f%%) | time %.2f s', nnz(mask_SPB_c), 100*nnz(mask_SPB_c)/numel(mask_SPB_c), toc(t));

% ISPB
t = tic;
for i = 1:numel(S_ispb)
    in = inpolygon(Xc, Yc, S_ispb(i).X/1000, S_ispb(i).Y/1000);
    mask_ISPB_c = mask_ISPB_c | in;
    if mod(i,50)==0 || i==numel(S_ispb), logf('ISPB %d/%d parts...', i, numel(S_ispb)); end
end
logf('ISPB mask crop nnz=%d (%.3f%%) | time %.2f s', nnz(mask_ISPB_c), 100*nnz(mask_ISPB_c)/numel(mask_ISPB_c), toc(t));

% OSPB
t = tic;
for i = 1:numel(S_ospb)
    in = inpolygon(Xc, Yc, S_ospb(i).X/1000, S_ospb(i).Y/1000);
    mask_OSPB_c = mask_OSPB_c | in;
    if mod(i,50)==0 || i==numel(S_ospb), logf('OSPB %d/%d parts...', i, numel(S_ospb)); end
end
logf('OSPB mask crop nnz=%d (%.3f%%) | time %.2f s', nnz(mask_OSPB_c), 100*nnz(mask_OSPB_c)/numel(mask_OSPB_c), toc(t));

%% ------------------ Stitch to full size (for QA / optional save) ----------
mask_SPB  = false(size(X_km));  mask_SPB( rmin:rmax, cmin:cmax) = mask_SPB_c;
mask_ISPB = false(size(X_km));  mask_ISPB(rmin:rmax, cmin:cmax) = mask_ISPB_c;
mask_OSPB = false(size(X_km));  mask_OSPB(rmin:rmax, cmin:cmax) = mask_OSPB_c;

logf('Full-size nnz: SPB=%d ISPB=%d OSPB=%d', nnz(mask_SPB), nnz(mask_ISPB), nnz(mask_OSPB));

%% ------------------ Save (cropped or full) --------------------------
outdir  = 'datasets_for_gmin';
if ~exist(outdir,'dir'), mkdir(outdir); end
outfile = fullfile(outdir,'spb_masks.mat');

t = tic;
if useCroppedSave
    % Save cropped masks + indices (smallest). Use sparse if very few trues.
    tosparse = @(A) (nnz(A) < sparseThreshold*numel(A));
    if tosparse(mask_SPB_c),  mask_SPB_c  = sparse(mask_SPB_c);  end
    if tosparse(mask_ISPB_c), mask_ISPB_c = sparse(mask_ISPB_c); end
    if tosparse(mask_OSPB_c), mask_OSPB_c = sparse(mask_OSPB_c); end
    % rect_mask is usually dense — keep dense for simplicity.
    save(outfile, 'mask_SPB_c','mask_ISPB_c','mask_OSPB_c', ...
                  'rect_mask','rmin','rmax','cmin','cmax','-v7.3');
    logf('Saved CROPPED masks → %s (%.2f s).', outfile, toc(t));
    logf('To reconstruct full masks later, stitch with rmin/rmax/cmin/cmax.');
else
    %save(outfile, 'mask_SPB','mask_ISPB','mask_OSPB','rect_mask','-v7');
    %logf('Saved FULL masks → %s (%.2f s).', outfile, toc(t));
end

%% ------------------ Quick QA figure ------------------------------
try
    t = tic;
    fig = figure('Visible','off','Color','w','Position',[100 100 900 800]);
    ax = axes('Parent',fig); hold(ax,'on'); axis(ax,'image'); box(ax,'on'); set(ax,'YDir','normal');
    imagesc(ax, X_km(1,:), Y_km(:,1), rect_mask); colormap(ax, gray);
    contour(ax, X_km(1,:), Y_km(:,1), mask_ISPB, [0.5 0.5], 'r', 'LineWidth',1.2);
    contour(ax, X_km(1,:), Y_km(:,1), mask_OSPB, [0.5 0.5], 'b', 'LineWidth',1.2);
    legend(ax, {'Rect','ISPB','OSPB'}, 'Location','southoutside'); xlabel(ax,'X (km)'); ylabel(ax,'Y (km)');
    qa_path = fullfile(outdir,'qa_spb_masks.png');
    exportgraphics(fig, qa_path, 'Resolution', 220);
    close(fig);
    logf('QA figure saved: %s (%.2f s)', qa_path, toc(t));
catch ME
    logf('QA figure failed: %s', ME.message);
end

logf('Done in %.2f s total.', toc(t0));
end

% ================= Helpers =================
function M = build_rect_mask(X_km, Y_km, rectCfg)
    M = true(size(X_km));  % default all pass
    if nargin<3 || ~isstruct(rectCfg), return; end
    if ~isfield(rectCfg,'enable') || ~rectCfg.enable, return; end
    if isfield(rectCfg,'rects') && ~isempty(rectCfg.rects)
        rects = rectCfg.rects;
    elseif all(isfield(rectCfg,{'xlim','ylim'})) && ~isempty(rectCfg.xlim) && ~isempty(rectCfg.ylim)
        rects = [rectCfg.xlim(:).' rectCfg.ylim(:).'];
    else
        return
    end
    if size(rects,2) ~= 4, error('cfg.rect.rects must be N×4: [xmin xmax ymin ymax]'); end
    keep = false(size(X_km));
    for r = 1:size(rects,1)
        xmin = min(rects(r,1), rects(r,2)); xmax = max(rects(r,1), rects(r,2));
        ymin = min(rects(r,3), rects(r,4)); ymax = max(rects(r,3), rects(r,4));
        keep = keep | (X_km >= xmin & X_km <= xmax & Y_km >= ymin & Y_km <= ymax);
    end
    mode = 'include';
    if isfield(rectCfg,'mode') && ~isempty(rectCfg.mode), mode = lower(string(rectCfg.mode)); end
    if strcmp(mode,'exclude'), M = ~keep; else, M = keep; end
end
