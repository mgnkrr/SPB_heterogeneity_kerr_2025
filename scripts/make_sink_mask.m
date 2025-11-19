function make_sink_mask()
% MAKE_SINK_MASK
% Build sink mask from hydraulic head, then keep only sink pixels that
% pass the *valid specularity* filter (Q ~= 9999 & finite) BEFORE forming
% connected components. Saves both pixel mask and pruned component labels.

%% -------------------- CONFIG --------------------
conn        = 8;                      % 4 or 8 connectivity
min_size    = 5;                      % min pixels per component
out_sinks   = '/Users/megankerr/Documents/Location-oldest-ice/datasets_for_gmin/sink_mask_new.mat';
out_comp    = '/Users/megankerr/Documents/Location-oldest-ice/datasets_for_gmin/sink_mask_comp.mat';
addpath('/Users/megankerr/Documents/Location-oldest-ice/data_ipics/Bed Head v1.0/Bed Head');

%% -------------------- Load geoid-corrected rasters --------------------
disp('Loading input data...');
load('datasets_for_gmin/coldex_bedelv.mat', 'B', 'Xgrid', 'Ygrid');      % B, Xgrid, Ygrid
load('datasets_for_gmin/REMA_srfelv.mat', 'REMA_srfelv');                % optional alt surface
load('datasets_for_gmin/coldex_srfelv.mat', 'S');                        % S (surface)

assert(isequal(size(B),size(REMA_srfelv),size(S),size(Xgrid),size(Ygrid)), ...
    'Size mismatch among rasters');

% Basic diagnostics on grid orientation & spacing
xv = Xgrid(1,:);  yv = Ygrid(:,1);
dx = median(diff(xv)); dy = median(diff(yv));
fprintf('[grid] size = %dx%d | X %s (dx=%.3f) | Y %s (dy=%.3f)\n', ...
    size(B,1), size(B,2), tern(diff(xv)>0,'asc','desc'), dx, tern(diff(yv)>0,'asc','desc'), dy);
assert(all(abs(diff(xv) - dx) < 10*eps(max(abs(xv)))),'X spacing not uniform');
assert(all(abs(diff(yv) - dy) < 10*eps(max(abs(yv)))),'Y spacing not uniform');

%% -------------------- Ice thickness --------------------
proj = projcrs(3031); 
[lat, ~] = projinv(proj, Xgrid, Ygrid); 
% lon unused % 
mask_coldex = (lat <= -88); 
% use COLDEX near the pole 
combined_srfelv = REMA_srfelv; 
combined_srfelv(mask_coldex) = S(mask_coldex); 
% % % Ice thickness (nonnegative) 
thk = max(combined_srfelv - B, 0);

%% -------------------- Normalize to north-up for GRIDobj ----------------
% GRIDobj expects y to be strictly decreasing (north -> south) down rows.
flipY = yv(1) < yv(end);  % true if south-up (ascending y)
if flipY
    disp('[orient] Detected south-up input (Y ascending). Flipping to north-up for processing.');
    Bz   = flipud(B);
    Thkz = flipud(thk);
    yv_north = flipud(yv);     % now strictly decreasing
else
    Bz   = B;
    Thkz = thk;
    yv_north = yv;             % already decreasing
end
xv_east = xv;                  % assume X increases eastward (common in EPSG:3031)

%% -------------------- GRIDobjs for TopoToolbox ------------------------
disp('Creating GRID objects (north-up) ...');
bed_dem  = GRIDobj(xv_east, yv_north, Bz);
thk_dem  = GRIDobj(xv_east, yv_north, Thkz);

%% -------------------- Hydraulic head --------------------
try
    disp('Calculating hydraulic head with bedhead() ...');
    head_dem = bedhead('bed', bed_dem, 'thickness', thk_dem);
catch
    warning('bedhead() not found; using bed elevation as head proxy');
    head_dem = bed_dem;
end

%% -------------------- Fill sinks & base sink pixels -------------------
disp('Computing DEM fill (north-up) ...');
DEMf = fillsinks(head_dem);

sinks_depth_north = DEMf.Z - head_dem.Z;           % north-up orientation
sink_mask_north   = (sinks_depth_north > 0) & isfinite(sinks_depth_north);

% Flip back to original orientation for saving/plotting
if flipY
    sink_mask = flipud(sink_mask_north);
else
    sink_mask = sink_mask_north;
end

%% -------------------- Save sink pixel mask ----------------------------
if ~exist(fileparts(out_sinks), 'dir'); mkdir(fileparts(out_sinks)); end
save(out_sinks, 'sink_mask', 'Xgrid', 'Ygrid', '-v7.3');
fprintf('Saved sink mask to: %s\n', out_sinks);

%% -------------------- Connected components (after filter) -------------
disp('Finding connected components ...');
CC = bwconncomp(sink_mask, conn);
fprintf('  initial components (post-filter): %d\n', CC.NumObjects);

keepC = cellfun(@numel, CC.PixelIdxList) >= min_size;
sink_mask_comp = false(size(sink_mask));
if any(keepC)
    sink_mask_comp(vertcat(CC.PixelIdxList{keepC})) = true;
end

CC2 = bwconncomp(sink_mask_comp, conn);
comp_id = uint32(labelmatrix(CC2));
fprintf('  kept components: %d (removed %d)\n', CC2.NumObjects, CC.NumObjects - CC2.NumObjects);

%% -------------------- Save component labels ---------------------------
if ~exist(fileparts(out_comp), 'dir'); mkdir(fileparts(out_comp)); end
save(out_comp, 'sink_mask_comp', 'comp_id', 'Xgrid', 'Ygrid', '-v7');
fprintf('Saved component mask/labels -> %s\n', out_comp);

%% -------------------- Quicklook (components) --------------------------
figure('Color','w','Name','Kept components over bed');
imagesc(xv, yv, B); axis image; set(gca,'YDir','normal'); colormap gray; hold on
title(sprintf('Kept sink components (n=%d)', CC2.NumObjects));
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on
contour(Xgrid, Ygrid, double(sink_mask_comp), [0.5 0.5], 'm', 'LineWidth', 2);

end

% --- tiny helper ---
function s = tern(cond, a, b)
if all(cond), s = a; else, s = b; end
end
