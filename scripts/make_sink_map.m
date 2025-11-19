function make_sink_map(S, varargin)
% MAKE_SINK_MAP  Professional sink-component map for publication.
%
% Features
% --------
% • Grayscale ice-thickness base (muted, annotation-friendly)
% • Optional residual overlay (separate axes + alpha by validity)
% • Sink components: wet = white fill (truecolor alpha), dry = none, black outlines (contour)
% • ISPB / OSPB outlines (meters shapefiles -> km plotting)
% • Inset scalebar (auto length selection)
% • Optional stats table (helper included; commented near bottom)
% • High-DPI export (PNG) and optional vector/bitmap PDF
% • Shapefile export of sink components with {ID,isWet,nPix,Area_km2}
%
% Example
% -------
%   L = load('figs_out/artifacts/S_plot_YYYYMMDD_HHMMSS.mat');  % S_plot.*
%   make_sink_map(L.S_plot, ...
%       'outdir','maps_out', 'run_id','rerun', ...
%       'ispb_shp','datasets_for_gmin/inner_spb.shp', ...
%       'ospb_shp','datasets_for_gmin/outer_spb.shp', ...
%       'residual',true, 'visible',true, 'export_pdf',true);

% ---------------- Options ----------------
p = inputParser;
addParameter(p,'outdir','maps_out',@(s)ischar(s)||isstring(s));
addParameter(p,'run_id','',@(s)ischar(s)||isstring(s));
addParameter(p,'visible',false,@(x)islogical(x)||ismember(x,[0 1]));
addParameter(p,'fontsize',12,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'fontname','Times New Roman',@(s)ischar(s)||isstring(s));
addParameter(p,'residual',true,@islogical);               % overlay residual if present
addParameter(p,'contours',[1000 1500 2000 2500 3000 3500 4000],@(v)isnumeric(v));
addParameter(p,'dpi',300,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'export_pdf',false,@islogical);
addParameter(p,'pdf_vector',true,@islogical);             % false -> bitmap PDF

% Region shapefiles (meters)
addParameter(p,'ispb_shp','',@(s)ischar(s)||isstring(s));
addParameter(p,'ospb_shp','',@(s)ischar(s)||isstring(s));

parse(p,varargin{:});
opt = p.Results;
if ~exist(opt.outdir,'dir'), mkdir(opt.outdir); end

% ---------------- Required fields ----------------
req = {'H','X_km','Y_km','sink_mask_comp','roi_mask','y_raw_full'};
for k = 1:numel(req)
    assert(isfield(S,req{k}), 'S.%s is required', req{k});
end
assert(isequal(size(S.H), size(S.X_km), size(S.Y_km)), 'Grid/H size mismatch');
assert(isequal(size(S.sink_mask_comp), size(S.H)), 'sink_mask_comp vs H size mismatch');
assert(isequal(size(S.roi_mask),       size(S.H)), 'roi_mask vs H size mismatch');
assert(isequal(size(S.y_raw_full),     size(S.H)), 'y_raw_full vs H size mismatch');

has_resid = isfield(S,'resid_db') && ~isempty(S.resid_db) && all(size(S.resid_db)==size(S.H));

% ---------------- Sink components (ROI-scoped) ----------------
mask_full = logical(S.roi_mask) & logical(S.sink_mask_comp);
if ~any(mask_full(:))
    warning('No sink components within ROI. Background-only map will be produced.');
end
CC = bwconncomp(mask_full, 8);

% Component wet/dry classification (any wet pixel -> wet component)
Yfull = logical(S.y_raw_full);
wet_comp = false(CC.NumObjects,1);
for i = 1:CC.NumObjects
    wet_comp(i) = any(Yfull(CC.PixelIdxList{i}));
end

% Build a wet-only mask (world-size)
wet_mask = false(size(mask_full));
for i = 1:CC.NumObjects
    if wet_comp(i)
        wet_mask(CC.PixelIdxList{i}) = true;
    end
end

% ---------------- Figure & base ----------------
fig = figure('Color','w','Name','Sink Components Map', ...
             'Visible', tern(opt.visible,'on','off'), ...
             'Units','pixels','Position',[80 80 1250 930]);

axH = axes('Parent',fig);
set(axH, 'FontName',opt.fontname, 'FontSize',opt.fontsize, ...
         'LineWidth',0.75, 'TickDir','out', 'Box','on', ...
         'XAxisLocation','bottom','YAxisLocation','left', ...
         'XColor','k','YColor','k');

% Hide toolbar in exports
try, axtoolbar(axH,'none'); catch, end

% Muted grayscale base
n = 256; lo = 0.96; hi = 0.72;
cmH = lo + (hi-lo)*gray(n);
imagesc(axH, S.X_km(1,:), S.Y_km(:,1), S.H);
axis(axH,'image'); set(axH,'YDir','normal'); hold(axH,'on');
colormap(axH, cmH);
cbH = colorbar(axH,'eastoutside'); cbH.Label.String = 'ice thickness (m)';
cbH.FontSize = max(8, opt.fontsize-1);

% Subtle thickness contours
if ~isempty(opt.contours)
    contour(axH, S.X_km, S.Y_km, S.H, opt.contours, ...
        'LineColor',[0.42 0.42 0.42], 'LineWidth',0.35, 'ShowText','off');
end

% ---------------- Residual overlay (separate axes) ----------------
if opt.residual && has_resid
    axR = axes('Parent',fig,'Position',get(axH,'Position'));
    try, axtoolbar(axR,'none'); catch, end
    hR  = imagesc(axR, S.X_km(1,:), S.Y_km(:,1), S.resid_db);
    set(axR,'YDir','normal','Color','none','Box','off'); axis(axR,'image');
    A = 0.55 * isfinite(S.resid_db); set(hR,'AlphaData', A);
    colormap(axR, parula);
    rlo = prctile(S.resid_db(isfinite(S.resid_db)), 2);
    rhi = prctile(S.resid_db(isfinite(S.resid_db)),98);
    if isfinite(rlo) && isfinite(rhi) && rhi>rlo, caxis(axR,[rlo rhi]); end
    cbR = colorbar(axR,'eastoutside'); cbR.Label.String = 'residual (dB)';
    cbR.FontSize = max(8,opt.fontsize-1);
    set(axR,'XAxisLocation','top','YAxisLocation','right','Visible','off','XTick',[],'YTick',[]);
    linkaxes([axH axR],'xy');
    axes(axR); hold(axR,'on');  % draw polygons above overlay
else
    axR = axH;
end

% ---------------- ISPB / OSPB outlines ----------------
hISPB = []; hOSPB = [];
if ~isempty(opt.ospb_shp) && isfile(opt.ospb_shp)
    S_ospb = shaperead(opt.ospb_shp);
    for k = 1:numel(S_ospb)
        plot(axR, S_ospb(k).X/1000, S_ospb(k).Y/1000, ...
            'LineStyle','-', 'LineWidth', 2.0, 'Color','k', 'HandleVisibility','off');
    end
    hOSPB = plot(axR, NaN,NaN,'-','LineWidth',2.0,'Color','k','DisplayName','Outer SPB');
end
if ~isempty(opt.ispb_shp) && isfile(opt.ispb_shp)
    S_ispb = shaperead(opt.ispb_shp);
    for k = 1:numel(S_ispb)
        plot(axR, S_ispb(k).X/1000, S_ispb(k).Y/1000, ...
            'LineStyle','--','LineWidth', 2.0, 'Color','k', 'HandleVisibility','off');
    end
    hISPB = plot(axR, NaN,NaN,'--','LineWidth',2.0,'Color','k','DisplayName','Inner SPB');
end

% ---------------- Sink polygons (world coords; orientation-safe) ---------
% 1) Fill wet components with a white truecolor image masked by AlphaData.
whiteRGB = ones([size(wet_mask), 3]);  % truecolor white, independent of colormap
hWetImg = image(axR, 'XData', S.X_km(1,:), 'YData', S.Y_km(:,1), 'CData', whiteRGB);
set(axR,'YDir','normal'); axis(axR,'image');   % ensure consistent orientation
set(hWetImg, 'AlphaData', double(wet_mask));   % 1 inside wet comps, 0 elsewhere
uistack(hWetImg, 'top');

% 2) Draw black outlines for all components using contour on the binary mask.
hold(axR,'on');
contour(axR, S.X_km, S.Y_km, double(mask_full), [0.5 0.5], ...
        'LineColor','k', 'LineWidth', 1.0);

% Legend proxies (appearance matches fill/edge)
hWet = patch('Parent',axR, 'XData',nan,'YData',nan, ...
             'FaceColor','w','FaceAlpha',1.0, ...
             'EdgeColor','k','LineWidth',1.0, ...
             'DisplayName','Wet sink');
hDry = patch('Parent',axR, 'XData',nan,'YData',nan, ...
             'FaceColor','none','FaceAlpha',1.0, ...
             'EdgeColor','k','LineWidth',1.0, ...
             'DisplayName','Dry sink');

% ---------------- Axis cosmetics ----------------
adjust_ticks(axR);
xlabel(axR,'Easting (km)','FontSize',opt.fontsize,'FontName',opt.fontname);
ylabel(axR,'Northing (km)','FontSize',opt.fontsize,'FontName',opt.fontname);
set(axR,'FontName',opt.fontname,'FontSize',opt.fontsize);
ylim(axR,[-200 250]);   % your desired window
xlim(axR,[0 700]);
grid(axR,'off');

% Legend (include SPB outlines if present)
legItems = [hWet hDry];
if ~isempty(hISPB), legItems = [legItems hISPB]; end
if ~isempty(hOSPB), legItems = [legItems hOSPB]; end
lgd = legend(axR, legItems, 'Orientation','horizontal', ...
    'NumColumns', numel(legItems), ...
    'Box','on', 'EdgeColor',[0.75 0.75 0.75], 'Color','w', 'Location','southoutside');
lgd.FontName = opt.fontname; lgd.FontSize = max(8,opt.fontsize-1);
lgd.ItemTokenSize = [18 9];

% ---------------- (Optional) Region stats table ----------------
%{
dx = double(abs(median(diff(S.X_km(1,:)), 'omitnan')));
dy = double(abs(median(diff(S.Y_km(:,1)), 'omitnan')));
pix_area_km2 = dx * dy;
polyMask = @(shpfile) build_poly_mask_km(shpfile, S.X_km, S.Y_km);
M_all  = true(size(S.H));
M_ispb = M_all; if ~isempty(opt.ispb_shp) && isfile(opt.ispb_shp), M_ispb = polyMask(opt.ispb_shp); end
M_ospb = M_all; if ~isempty(opt.ospb_shp) && isfile(opt.ospb_shp), M_ospb = polyMask(opt.ospb_shp); end
M_sink = mask_full;
stats_all  = region_stats(M_sink & M_all,  Yfull, pix_area_km2);
stats_ispb = region_stats(M_sink & M_ispb, Yfull, pix_area_km2);
stats_ospb = region_stats(M_sink & M_ospb, Yfull, pix_area_km2);
headers = {'Region','N Pixels','N Comps','Avg Area (km^2)','Wet %','Dry %'};
rows = {
  'ISPB', commafy(stats_ispb.n_pixels), commafy(stats_ispb.n_comps), numfmt(stats_ispb.avg_area_km2), numfmt(stats_ispb.pct_wet), numfmt(stats_ispb.pct_dry);
  'OSPB', commafy(stats_ospb.n_pixels), commafy(stats_ospb.n_comps), numfmt(stats_ospb.avg_area_km2), numfmt(stats_ospb.pct_wet), numfmt(stats_ospb.pct_dry);
  'All',  commafy(stats_all.n_pixels),  commafy(stats_all.n_comps),  numfmt(stats_all.avg_area_km2),  numfmt(stats_all.pct_wet),  numfmt(stats_all.pct_dry);
};
posMap     = get(axR,'Position');
tableWidth = posMap(3) * 0.60;
tableHeight= posMap(4) * 0.13;
tableY     = posMap(2) + posMap(4) - 0.064;
tableX     = posMap(1) + (posMap(3)-tableWidth)/2;
axTbl = axes('Parent',fig,'Units','normalized', ...
             'Position',[tableX tableY tableWidth tableHeight], ...
             'Color','w','Box','on');
axis(axTbl,'off');
draw_stats_table_opaque(axTbl, headers, rows, max(8, opt.fontsize-1), opt.fontname, 'altrows',true);
uistack(axTbl,'top');
%}

% ---------------- Orientation sanity check (quiet) ----------------
try
    if has_resid
        R0  = double(S.resid_db); R1 = flipud(R0);
        M   = logical(S.roi_mask) & isfinite(R0) & isfinite(Yfull);
        thr = prctile(R0(M), 80);
        s0 = mean( (R0(M) >= thr) == Yfull(M) );
        s1 = mean( (R1(M) >= thr) == Yfull(M) );
        if s1 > s0 + 0.02
            warning('Residual appears improved when vertically flipped (check orientation).');
        end
    end
catch
end

% ---------------- Shapefile export: sink components (wet/dry) -------------
% One record per connected component inside ROI.
% Geometry: polygons (meters). Attributes: centroid + wet flags + counts.
try
    dx = double(abs(median(diff(S.X_km(1,:)), 'omitnan')));
    dy = double(abs(median(diff(S.Y_km(:,1)), 'omitnan')));
    pix_area_km2 = dx * dy;

    tmpl = struct('Geometry','Polygon', ...
                  'X',double([]), 'Y',double([]), ...
                  'ID',double(0), ...
                  'CentroidX',double(0), 'CentroidY',double(0), ...
                  'isWetNum',double(0), 'isWetStr','', ...
                  'nPix',double(0), 'Area_km2',double(0));
    recs = repmat(tmpl, CC.NumObjects, 1);

    for i = 1:CC.NumObjects
        Mi = false(size(mask_full));
        Mi(CC.PixelIdxList{i}) = true;

        Bi = bwboundaries(Mi, 8, 'noholes');

        Xall = []; Yall = [];
        for j = 1:numel(Bi)
            b  = Bi{j}; rr = b(:,1); cc = b(:,2);
            x = double(S.X_km(1,cc)) * 1000;  % km -> m
            y = double(S.Y_km(rr,1)) * 1000;
            if x(1)~=x(end) || y(1)~=y(end)
                x(end+1,1) = x(1);
                y(end+1,1) = y(1);
            end
            Xall = [Xall; x(:); NaN]; %#ok<AGROW>
            Yall = [Yall; y(:); NaN]; %#ok<AGROW>
        end

        Rprops = regionprops(Mi, 'Centroid');   % pixel coords [col, row]
        if ~isempty(Rprops)
            cxy_pix = Rprops(1).Centroid;
            cX_m = double(S.X_km(1, round(cxy_pix(1)))) * 1000;
            cY_m = double(S.Y_km(round(cxy_pix(2)), 1)) * 1000;
        else
            v = ~isnan(Xall) & ~isnan(Yall);
            cX_m = mean(Xall(v), 'omitnan');
            cY_m = mean(Yall(v), 'omitnan');
        end

        npix = double(numel(CC.PixelIdxList{i}));
        wet  = double(wet_comp(i));

        recs(i).X         = Xall;
        recs(i).Y         = Yall;
        recs(i).ID        = double(i);
        recs(i).CentroidX = cX_m;
        recs(i).CentroidY = cY_m;
        recs(i).isWetNum  = wet;                     % 0/1 for DBF
        recs(i).isWetStr  = tern(logical(wet),'WET','DRY');
        recs(i).nPix      = npix;
        recs(i).Area_km2  = npix * pix_area_km2;
    end

    shpfile = fullfile(opt.outdir, ['sinks_' char(opt.run_id) '.shp']);
    shapewrite(recs, shpfile);
    fprintf('[make_sink_map] Saved shapefile -> %s (records=%d)\n', shpfile, numel(recs));

catch
    warning('Shapefile export failed');
end

% ---------------- Export ----------------
runfrag = ''; if ~isempty(opt.run_id), runfrag = ['_' char(opt.run_id)]; end
pngfile = fullfile(opt.outdir, ['map_sinks' runfrag '.png']);
exportgraphics(fig, pngfile, 'Resolution', opt.dpi);
fprintf('[make_sink_map] Saved PNG -> %s | components: wet=%d, dry=%d\n', pngfile, nnz(wet_comp), CC.NumObjects-nnz(wet_comp));

if opt.export_pdf
    pdffile = fullfile(opt.outdir, ['map_sinks' runfrag '.pdf']);
    if opt.pdf_vector
        exportgraphics(fig, pdffile, 'ContentType','vector');
    else
        exportgraphics(fig, pdffile, 'ContentType','image', 'Resolution', opt.dpi);
    end
    fprintf('[make_sink_map] Saved PDF  -> %s\n', pdffile);
end

if ~opt.visible, close(fig); end
end

% ===================== Helpers =====================

function adjust_ticks(ax)
% Make ticks readable and not crowded; aim for ~6 ticks per axis.
    xl = xlim(ax); yl = ylim(ax);
    nx = nice_count(range(xl), 6);
    ny = nice_count(range(yl), 6);
    set(ax, 'XTickMode','auto', 'YTickMode','auto');
    xt = linspace(xl(1), xl(2), nx);
    yt = linspace(yl(1), yl(2), ny);
    xt = neat_round(xt); yt = neat_round(yt);
    set(ax,'XTick',xt,'YTick',yt);
end

function n = nice_count(span, target)
    if span<=0 || ~isfinite(span), n = target; return; end
    n = max(3, min(8, target));
end

function v = neat_round(v)
    if all(~isfinite(v)), return; end
    span = max(v)-min(v);
    if span==0 || ~isfinite(span), return; end
    step = span / max(1,(numel(v)-1));
    p10 = 10^floor(log10(step));
    stepN = step/p10;
    if stepN<1.5, s=1; elseif stepN<3, s=2; elseif stepN<7, s=5; else, s=10; end
    s = s*p10;
    v = round(v/s)*s;
end

function r = range(x), r = max(x)-min(x); end

function draw_scalebar(ax, varargin)
% Draw a simple scalebar in axes normalized coordinates (units: km).
    ip = inputParser;
    addParameter(ip,'anchor','sw',@(s)ischar(s)||isstring(s));
    addParameter(ip,'pad',[0.06 0.06],@(v)isnumeric(v)&&numel(v)==2);
    addParameter(ip,'height',0.015,@(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip,varargin{:});
    o = ip.Results;

    frac = 0.20;  % 20% of axis width
    xl = xlim(ax); yl = ylim(ax);
    Lkm_target = frac * (xl(2)-xl(1));
    L = nice_length(Lkm_target);
    Lnorm = L / (xl(2)-xl(1));

    switch lower(o.anchor)
        case 'sw', x0=0+o.pad(1); y0=0+o.pad(2);
        case 'se', x0=1-o.pad(1)-Lnorm; y0=0+o.pad(2);
        case 'nw', x0=0+o.pad(1); y0=1-o.pad(2)-o.height;
        case 'ne', x0=1-o.pad(1)-Lnorm; y0=1-o.pad(2)-o.height;
        otherwise, x0=0.06; y0=0.06;
    end

    x1 = xl(1) + x0 * (xl(2)-xl(1));
    x2 = x1 + L;
    y1 = yl(1) + y0 * (yl(2)-yl(1));
    h  = o.height * (yl(2)-yl(1));

    hold(ax,'on');
    rectangle(ax,'Position',[x1 y1 (x2-x1) h], 'FaceColor','k', 'EdgeColor','k', 'LineWidth',0.75);
    rectangle(ax,'Position',[x1 y1+h*0.40 (x2-x1) h*0.20], 'FaceColor','w', 'EdgeColor','none');

    txt = sprintf('%g km', L);
    text(ax, (x1+x2)/2, y1 - 0.8*h, txt, 'HorizontalAlignment','center', ...
         'VerticalAlignment','top', 'FontName', get(ax,'FontName'), ...
         'FontSize', max(8, get(ax,'FontSize')-1), 'Color','k');
end

function L = nice_length(Lt)
    if Lt<=0 || ~isfinite(Lt), L=1; return; end
    p10 = 10^floor(log10(Lt));
    n = Lt/p10;
    if n<1.5, b=1; elseif n<3.5, b=2; elseif n<7.5, b=5; else, b=10; end
    L = b*p10;
end

function M = build_poly_mask_km(shpfile, X_km, Y_km)
    Sshp = shaperead(shpfile);  % X/Y in meters
    M = false(size(X_km));
    for k = 1:numel(Sshp)
        xk = Sshp(k).X / 1000;  yk = Sshp(k).Y / 1000;
        isn = isnan(xk) | isnan(yk);
        if any(isn)
            ii = [0 find(isn) numel(xk)+1];
            for j = 1:numel(ii)-1
                seg = (ii(j)+1):(ii(j+1)-1);
                if isempty(seg), continue; end
                M = M | inpolygon(X_km, Y_km, xk(seg), yk(seg));
            end
        else
            M = M | inpolygon(X_km, Y_km, xk, yk);
        end
    end
end

function st = region_stats(M, Yfull, pix_area_km2)
    M = M & isfinite(Yfull);
    n_pixels = nnz(M);
    CC = bwconncomp(M, 8);
    n_comps = CC.NumObjects;
    if n_comps > 0
        comp_sizes = cellfun(@numel, CC.PixelIdxList);
        comp_area_km2 = comp_sizes * pix_area_km2;
        avg_area = mean(comp_area_km2);
        wet_comp = false(n_comps,1);
        for i = 1:n_comps
            wet_comp(i) = any(Yfull(CC.PixelIdxList{i}));
        end
        pct_wet = 100 * nnz(wet_comp) / n_comps;
        pct_dry = 100 - pct_wet;
    else
        avg_area = NaN; pct_wet = NaN; pct_dry = NaN;
    end
    st = struct('n_pixels',n_pixels,'n_comps',n_comps, ...
                'avg_area_km2',avg_area,'pct_wet',pct_wet,'pct_dry',pct_dry);
end

function draw_stats_table_opaque(ax, headers, rows, basefs, fontname, varargin)
    ip = inputParser; addParameter(ip,'altrows',true,@islogical); parse(ip,varargin{:});
    alt = ip.Results.altrows;
    cla(ax); hold(ax,'on'); axis(ax,'off'); set(ax,'XLim',[0 1],'YLim',[0 1]);
    nr = size(rows,1) + 1;  nc = numel(headers);
    pad = 0.012; x0=0; y0=0; w=1; h=1;
    colw = w * ones(1,nc) / nc;  rowh = h * ones(1,nr) / nr;
    cGrid=[.80 .80 .80]; cHdr=[.94 .94 .94]; cEdge=[.75 .75 .75]; cAlt=[.985 .985 .985];

    rectangle('Parent',ax,'Position',[x0,1-rowh(1),w,rowh(1)], ...
              'FaceColor',cHdr,'EdgeColor',cEdge,'LineWidth',0.75);

    for r = 2:nr
        fc = 'w'; if alt && mod(r,2)==0, fc = cAlt; end
        rectangle('Parent',ax,'Position',[x0,1-sum(rowh(1:r)),w,rowh(r)], ...
                  'FaceColor',fc,'EdgeColor',cEdge,'LineWidth',0.5);
    end

    for c = 1:nc-1
        x = x0 + sum(colw(1:c));
        line(ax,[x x],[y0 y0+h],'Color',cGrid,'LineWidth',0.5);
    end
    for r = 1:nr-1
        y = 1 - sum(rowh(1:r));
        line(ax,[x0 x0+w],[y y],'Color',cGrid,'LineWidth',0.5);
    end

    fsH = max(8, basefs);
    for c = 1:nc
        cx = x0 + sum(colw(1:c)) - colw(c);
        cy = 1 - sum(rowh(1:1));
        tx = cx + pad; ty = cy + rowh(1)/2;
        text(ax, tx, ty, headers{c}, 'FontName',fontname,'FontWeight','bold', ...
             'FontSize', fsH, 'Color',[.05 .05 .05], 'HorizontalAlignment','left', ...
             'VerticalAlignment','middle', 'Clipping','on');
    end

    fsB = max(8, basefs-1);
    for r = 1:size(rows,1)
        for c = 1:nc
            val = rows{r,c};
            cx = x0 + sum(colw(1:c)) - colw(c);
            cy = 1 - sum(rowh(1:r+1));
            if c == 1, tx = cx + pad; ha='left'; else, tx = cx + colw(c) - pad; ha='right'; end
            ty = cy + rowh(r+1)/2;
            text(ax, tx, ty, val, 'FontName',fontname,'FontSize',fsB,'Color',[.1 .1 .1], ...
                 'HorizontalAlignment',ha,'VerticalAlignment','middle','Clipping','on');
        end
    end
end

function out = tern(cond,a,b), if cond, out=a; else, out=b; end, end
function s = commafy(x), if isnan(x), s='NaN'; else, s=regexprep(sprintf('%.0f',x),'(\d)(?=(\d{3})+(?!\d))','$1,'); end, end
function s = numfmt(x), if isnan(x), s='NaN'; elseif x>=100, s=sprintf('%.0f',x); elseif x>=10, s=sprintf('%.1f',x); else, s=sprintf('%.2f',x); end, end
