function cellImgs = segmentCells3Dbatch()
fld = uigetdir;
files = dir(fullfile(fld, '**', '*.ims'));
cellImgs = {};
for i=1:length(files)
    fname = fullfile(files(i).folder,files(i).name);

    dsetPath = sprintf('/DataSet/ResolutionLevel %d/TimePoint 0/Channel %d/Data', ...
        0, 0);

    % Read bit depth
    dataInfo  = h5info(fname, dsetPath);
    bitDepth  = dataInfo.Datatype.Size * 8;
    dataClass = bitDepthToClass(bitDepth);
    fprintf('Detected bit depth: %d-bit (%s)\n', bitDepth, dataClass);

    % Read volume
    fprintf('Reading volume...\n');
    tic;
    vol = h5read(fname, dsetPath);
    fprintf('Done. Read time: %.1f s | Size: %s | Class: %s\n', ...
        toc, mat2str(size(vol)), class(vol));

    % vol = bfOpen3DVolume(vol);
    % vol = vol{1}{1};
    out = segmentCells3D(vol);
    cellImgs = [cellImgs, out];
end
end
function cls = bitDepthToClass(bitDepth)
    switch bitDepth
        case 1,  cls = 'logical';
        case 8,  cls = 'uint8';
        case 16, cls = 'uint16';
        case 32, cls = 'uint32';
        case 64, cls = 'uint64';
        otherwise
            warning('Unrecognised bit depth %d, defaulting to uint16', bitDepth);
            cls = 'uint16';
    end
end