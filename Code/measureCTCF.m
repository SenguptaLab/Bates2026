function [CTCFs] = measureCTCF()

folder = uigetdir;
files  = dir(fullfile(folder, '*.ims'));
CTCFs  = [];

for f = 1:length(files)

    fname = fullfile(files(f).folder, files(f).name);

    dsetPath = sprintf('/DataSet/ResolutionLevel %d/TimePoint 0/Channel %d/Data', 0, 0);

    dataInfo  = h5info(fname, dsetPath);
    bitDepth  = dataInfo.Datatype.Size * 8;
    dataClass = bitDepthToClass(bitDepth);
    fprintf('Detected bit depth: %d-bit (%s)\n', bitDepth, dataClass);

    fprintf('Reading volume...\n');
    tic;
    img = h5read(fname, dsetPath);
    fprintf('Done. Read time: %.1f s | Size: %s | Class: %s\n', ...
            toc, mat2str(size(img)), class(img));

    img = double(img);
    if size(img,3) > 1
        img = max(img, [], 3);
    end

    row = process_one_file(img, files(f).name);
    CTCFs = [CTCFs; row];
end

end

% =========================================================================
function row = process_one_file(img, fname_display)

    vmin = prctile(img(:), 1);
    vmax = prctile(img(:), 99.9);

    vox          = img(:);
    vox_bright   = vox(vox > prctile(vox, 50));
    thresh_guess = mean(vox_bright) + std(vox_bright);
    thresh_guess = min(thresh_guess, double(vmax));

    % --- Threshold UI ---
    fig = uifigure('Name', fname_display);
    fig.Position = [50 50 1500 1000];

    ax = uiaxes(fig);
    ax.Position = [20 220 1460 750];
    imshow(img, [], 'Parent', ax, 'InitialMagnification', 'fit');
    hold(ax, 'on');

    sld = uislider(fig);
    sld.Limits   = [0 double(vmax)];
    sld.Value    = thresh_guess;
    sld.Position = [200 150 400 3];

    lbl_coarse = uilabel(fig);
    lbl_coarse.Position = [200 170 200 20];
    lbl_coarse.Text     = 'Threshold (coarse)';

    spn = uispinner(fig);
    spn.Limits   = [0 double(vmax)];
    spn.Value    = thresh_guess;
    spn.Step     = max(1, round((vmax - vmin) / 1000));
    spn.Position = [650 140 150 30];

    lbl_fine = uilabel(fig);
    lbl_fine.Position = [650 170 100 20];
    lbl_fine.Text     = 'Fine adjust';

    val_lbl = uilabel(fig);
    val_lbl.Position = [820 140 220 30];
    val_lbl.Text     = sprintf('Value: %.1f', thresh_guess);

    lbl_range = uilabel(fig);
    lbl_range.Position = [1060 170 80 20];
    lbl_range.Text     = 'Slider max';

    range_spn = uispinner(fig);
    range_spn.Limits   = [0 double(max(img(:)))];
    range_spn.Value    = double(vmax);
    range_spn.Step     = max(1, round(double(max(img(:))) / 1000));
    range_spn.Position = [1060 140 150 30];
    range_spn.ValueChangedFcn = @on_range_change;

    okBtn = uibutton(fig, 'push');
    okBtn.Text     = 'OK';
    okBtn.Position = [1230 140 100 30];

    sld.ValueChangedFcn   = @(src,~) on_slider(src, ax, img, val_lbl, spn);
    spn.ValueChangedFcn   = @(src,~) on_spinner(src, ax, img, val_lbl, sld);
    okBtn.ButtonPushedFcn = @(~,~) uiresume(fig);

    update_threshold_display(sld.Value, ax, img);

    uiwait(fig);

    mask = guidata(fig);
    mask = mask.mask;
    mask = bwareaopen(mask, 100);
    close(fig);

    % --- Region stats ---
    stats = regionprops(mask, img, 'PixelIdxList', 'Area', 'MeanIntensity', 'Centroid', 'Image', 'BoundingBox');
    n_obj = numel(stats);

    if n_obj == 0
        warning('No objects found in %s', fname_display);
        row = [0 0];
        return
    end

    cc        = bwconncomp(mask);
    label_img = labelmatrix(cc);
    img_disp  = (img - min(img(:))) / (max(img(:)) - min(img(:)));
    cmap      = lines(n_obj);
    selected  = false(1, n_obj);

    boundaries = cell(1, n_obj);
    for j = 1:n_obj
        obj_2d        = label_img == j;
        boundaries{j} = bwboundaries(obj_2d, 'noholes');
    end

    % --- Click-to-select figure ---
    fig2 = figure('Name', sprintf('%s — Click objects to select | Enter to confirm', fname_display), ...
                  'Color', 'k', 'KeyPressFcn', @on_key, ...
                  'WindowButtonDownFcn', @on_click);
    fig2.Position = [50 50 1400 1000];

    axsel = axes('Parent', fig2, 'Color', 'k');
    axis(axsel, 'image', 'off');

    draw_scene();

    uicontrol('Parent', fig2, 'Style', 'pushbutton', 'String', 'Confirm selection', ...
              'Units', 'normalized', 'Position', [0.35 0.01 0.3 0.05], ...
              'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'w', ...
              'FontSize', 11, 'Callback', @(~,~) set(fig2, 'UserData', 'done'));

    sel_lbl = uicontrol('Parent', fig2, 'Style', 'text', ...
                        'Units', 'normalized', 'Position', [0.01 0.01 0.3 0.04], ...
                        'String', 'Selected: none', ...
                        'BackgroundColor', 'k', 'ForegroundColor', 'w', ...
                        'FontSize', 10, 'HorizontalAlignment', 'left');

    waitfor(fig2, 'UserData', 'done');

    sel = find(selected);
    if ishandle(fig2); close(fig2); end

    if isempty(sel)
        row = [0 0];
        return
    end

    % --- Background ROI ---
    imSel  = figure('Name', 'Drag ROI to background then press Enter', ...
                    'Position', [50 50 1400 1000]);
    ax_roi = axes('Parent', imSel);
    imagesc(ax_roi, imadjust(mat2gray(img)));
    colormap(ax_roi, gray);
    axis(ax_roi, 'image', 'off');
    set(ax_roi, 'XDir', 'normal', 'YDir', 'reverse');
    hold(ax_roi, 'on');
    title(ax_roi, 'Drag ROI to background then press Enter', 'FontSize', 18);

    B   = boundaries{sel(1)};
    roi = images.roi.Polygon(ax_roi, 'Color', 'r');
    roi.Position            = [B{1}(:,2) B{1}(:,1)];
    roi.InteractionsAllowed = 'translate';
    pause;

    BKGmask = createMask(roi);
    BKG     = mean(img(BKGmask > 0), 'all');
    delete(roi);
    close(imSel);

    % --- CTCF calculation ---
    stats_sel = stats(sel);
    arr = [cat(1, stats_sel.MeanIntensity), cat(1, stats_sel.Area), repmat(BKG, numel(sel), 1)];
    row = [arr(:,1) .* arr(:,2) - (arr(:,3) .* arr(:,2)), arr(:,1)];

    % -------------------------------------------------------------------------
    function draw_scene()
        cla(axsel);
        imagesc('Parent', axsel, 'CData', img_disp);
        colormap(axsel, gray);
        axis(axsel, 'image', 'off');
        set(axsel, 'XDir', 'normal', 'YDir', 'reverse');
        hold(axsel, 'on');

        for jj = 1:n_obj
            col = cmap(jj,:);

            if selected(jj)
                obj_2d  = label_img == jj;
                overlay = zeros([size(img_disp) 3]);
                for c = 1:3
                    overlay(:,:,c) = col(c) * obj_2d;
                end
                image(axsel, 'CData', overlay, 'AlphaData', 0.55 * obj_2d);
                set(axsel, 'XDir', 'normal', 'YDir', 'reverse');
            end

            for b = 1:numel(boundaries{jj})
                bnd = boundaries{jj}{b};
                plot(axsel, bnd(:,2), bnd(:,1), '-', 'Color', col, 'LineWidth', 2);
            end

            cx = stats(jj).Centroid(1);
            cy = stats(jj).Centroid(2);

            % Mean intensity — shadow then white
            text(axsel, cx+1, cy + 13, sprintf('%.0f', stats(jj).MeanIntensity), ...
                 'Color', 'k', 'FontSize', 12, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
            text(axsel, cx, cy + 12, sprintf('%.0f', stats(jj).MeanIntensity), ...
                 'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end

        hold(axsel, 'off');
        title(axsel, 'Click objects to select  |  index / mean intensity  |  Enter to confirm', 'Color', 'w');
    end

    function on_click(~, ~)
        if gca ~= axsel; return; end
        pt     = get(axsel, 'CurrentPoint');
        col_pt = round(pt(1,1));
        row_pt = round(pt(1,2));

        col_pt = max(1, min(size(label_img,2), col_pt));
        row_pt = max(1, min(size(label_img,1), row_pt));

        jj = label_img(row_pt, col_pt);
        if jj == 0; return; end

        selected(jj) = ~selected(jj);
        draw_scene();

        sel_now = find(selected);
        if isempty(sel_now)
            sel_lbl.String = 'Selected: none';
        else
            sel_lbl.String = sprintf('Selected: %s', num2str(sel_now));
        end
    end

    function on_key(~, evt)
        if strcmp(evt.Key, 'return')
            set(fig2, 'UserData', 'done');
        end
    end

    function on_range_change(src, ~)
        new_max    = src.Value;
        sld.Limits = [0 new_max];
        spn.Limits = [0 new_max];
        if sld.Value > new_max
            sld.Value = new_max;
            spn.Value = new_max;
            val_lbl.Text = sprintf('Value: %.1f', new_max);
            update_threshold_display(new_max, ax, img);
        end
    end

end

% =========================================================================
function on_slider(src, ax, img, val_lbl, spn)
    t = src.Value;
    spn.Value    = t;
    val_lbl.Text = sprintf('Value: %.1f', t);
    update_threshold_display(t, ax, img);
end

% =========================================================================
function on_spinner(src, ax, img, val_lbl, sld)
    t = src.Value;
    sld.Value    = clamp(t, sld.Limits(1), sld.Limits(2));
    val_lbl.Text = sprintf('Value: %.1f', t);
    update_threshold_display(t, ax, img);
end

% =========================================================================
function update_threshold_display(t, ax, img)
    bw = img > t;
    bw = bwareaopen(bw, 100);

    cla(ax);
    imshow(img, [], 'Parent', ax, 'InitialMagnification', 'fit');
    hold(ax, 'on');
    B = bwboundaries(bw);
    for k = 1:length(B)
        boundary = B{k};
        plot(ax, boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5);
    end
    hold(ax, 'off');

    data = guidata(ax.Parent);
    data.mask = bw;
    guidata(ax.Parent, data);
end

% =========================================================================
function v = clamp(v, lo, hi)
    v = max(lo, min(hi, v));
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