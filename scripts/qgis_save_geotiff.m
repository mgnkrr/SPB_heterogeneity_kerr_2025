function qgis_save_geotiff(G, X, Y, outpath, nodata)
% QGIS_SAVE_GEOTIFF  Write G(x,y) to GeoTIFF with EPSG:3031 + NoData.
%   G, X, Y : full grids in meters (EPSG:3031), same size
%   outpath : e.g., 'maps_out/ghf_best.tif'
%   nodata  : e.g., -9999 (if omitted, NaNs remain as NoData)

    if nargin < 5, nodata = NaN; end
    [nrows, ncols] = size(G);
    assert(isequal(size(X),size(Y),[nrows,ncols]), 'X/Y size must match G');

    xmin = min(X(:)); xmax = max(X(:));
    ymin = min(Y(:)); ymax = max(Y(:));

    R = maprefcells([xmin xmax], [ymin ymax], [nrows ncols]);

    A = single(G);
    if ~isnan(nodata)
        A(~isfinite(A)) = nodata;
    end

    if ~endsWith(lower(outpath), {'.tif','.tiff'}), outpath = [outpath '.tif']; end

    ttags = struct();
    try
        ttags.Compression = Tiff.Compression.LZW;
        ttags.TileWidth   = 256;
        ttags.TileLength  = 256;
    catch
        ttags.Compression = 'LZW';
    end

    geotiffwrite(outpath, A, R, 'CoordRefSysCode', 3031, 'TiffTags', ttags);
    fprintf('[qgis] GeoTIFF written: %s (EPSG:3031)\n', outpath);
end
