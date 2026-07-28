function [cellImgs, stats] = segmentCells3D(vol)

fig = uifigure('Name', '3D Threshold Viewer');
fig.Position = [100 100 1076 717];

viewer = volshow(vol, 'Parent', fig);

vmin = prctile(vol(:), 1);
vmax = prctile(vol(:), 99.9);

% --- Best guess: mean + 1 std of bright voxels ---
vox       = double(vol(:));
vox_bright = vox(vox > prctile(vox, 50));  % upper half only
thresh_otsu = mean(vox_bright) + std(vox_bright);
thresh_otsu = min(thresh_otsu, double(vmax));  % clamp to slider range

% --- Coarse slider ---
lbl_coarse = uilabel(fig);
lbl_coarse.Position = [50 120 80 20];
lbl_coarse.Text     = 'Coarse';

sld = uislider(fig);
sld.Limits   = [0 800];
sld.Value    = thresh_otsu;
sld.Position = [130 120 250 3];

% --- Fine spinner ---
lbl_fine = uilabel(fig);
lbl_fine.Position = [390 105 40 20];
lbl_fine.Text     = 'Fine';

spn = uispinner(fig);
spn.Limits   = [0 800];
spn.Value    = thresh_otsu;
spn.Step     = max(1, round((vmax - vmin) / 1000));
spn.Position = [430 105 120 30];

% --- Value label ---
val_lbl = uilabel(fig);
val_lbl.Position = [560 105 160 30];
val_lbl.Text     = sprintf('Value: %.1f', thresh_otsu);

% --- Slider max spinner ---
lbl_range = uilabel(fig);
lbl_range.Position = [730 105 80 30];
lbl_range.Text     = 'Slider max';

range_spn = uispinner(fig);
range_spn.Limits   = [0 double(max(vol(:)))];
range_spn.Value    = double(vmax);
range_spn.Step     = max(1, round(double(max(vol(:))) / 1000));
range_spn.Position = [820 105 150 30];
range_spn.ValueChangedFcn = @on_range_change;

% --- OK button ---
okBtn = uibutton(fig, 'push');
okBtn.Text     = 'OK';
okBtn.Position = [990 105 60 30];

% --- Callbacks ---
sld.ValueChangedFcn   = @on_slider;
spn.ValueChangedFcn   = @on_spinner;
okBtn.ButtonPushedFcn = @(~,~) uiresume(fig);

% Apply initial threshold
updateVolume(sld, thresh_otsu, viewer, vol);

uiwait(fig);

mask = guidata(fig);
mask = mask.mask;
mask = imfill(mask, 'holes');
se   = strel('disk', 5);
mask = imclose(mask, se);
close(fig);

stats = regionprops3(mask, vol, 'all');
n_obj = size(stats, 1);

if n_obj == 0
    warning('No objects found after thresholding.');
    cellImgs = {};
    return
end

cc        = bwconncomp(mask);
label_vol = labelmatrix(cc);

% --- 2D selection on MIP ---
mip      = max(double(vol), [], 3);
mip_disp = (mip - min(mip(:))) / (max(mip(:)) - min(mip(:)));

cmap      = lines(n_obj);
selected  = false(1, n_obj);
label_mip = max(label_vol, [], 3);

% Precompute boundaries
boundaries = cell(1, n_obj);
for j = 1:n_obj
    obj_2d        = label_mip == j;
    boundaries{j} = bwboundaries(obj_2d, 'noholes');
end

fig2 = figure('Name', 'Click objects to select | Enter to confirm', ...
              'Color', 'k', 'KeyPressFcn', @on_key, ...
              'WindowButtonDownFcn', @on_click);
fig2.Position = [100 100 1000 800];

ax = axes('Parent', fig2, 'Color', 'k');
axis(ax, 'image', 'off');

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
    warning('No objects selected.');
    cellImgs = {};
    return
end

fprintf('Selected objects: %s\n', num2str(sel));

cellImgs = {};
for i = 1:numel(sel)
    idx         = sel(i);
    cell_mask   = stats.Image{idx};
    cell_double = NaN(size(cell_mask));
    cell_double(cell_mask) = double(vol(stats.VoxelIdxList{idx}));
    cellImgs{end+1} = cell_double;
end

% =========================================================================
    function draw_scene()
        cla(ax);
        imagesc('Parent', ax, 'CData', mip_disp);
        colormap(ax, gray);
        axis(ax, 'image', 'off');
        hold(ax, 'on');

        for j = 1:n_obj
            col = cmap(j,:);

            if selected(j)
                obj_2d  = label_mip == j;
                overlay = zeros([size(mip_disp) 3]);
                for c = 1:3
                    overlay(:,:,c) = col(c) * obj_2d;
                end
                image(ax, 'CData', overlay, 'AlphaData', 0.4 * obj_2d);
            end

            for b = 1:numel(boundaries{j})
                bnd = boundaries{j}{b};
                plot(ax, bnd(:,2), bnd(:,1), '-', ...
                     'Color', col, 'LineWidth', 1.5);
            end

            cx = stats.Centroid(j,1);
            cy = stats.Centroid(j,2);
            text(ax, cx, cy, num2str(j), 'Color', col, ...
                 'FontSize', 8, 'FontWeight', 'bold', ...
                 'HorizontalAlignment', 'center');
        end

        hold(ax, 'off');
        title(ax, 'Click objects to select | Enter to confirm', 'Color', 'w');
    end

    function on_click(~, ~)
        if gca ~= ax; return; end
        pt  = get(ax, 'CurrentPoint');
        col = round(pt(1,1));
        row = round(pt(1,2));

        col = max(1, min(size(label_mip,2), col));
        row = max(1, min(size(label_mip,1), row));

        j = label_mip(row, col);
        if j == 0; return; end

        selected(j) = ~selected(j);
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

    function on_slider(src, ~)
        t            = src.Value;
        spn.Value    = t;
        val_lbl.Text = sprintf('Value: %.1f', t);
        updateVolume(src, t, viewer, vol);
    end

    function on_spinner(src, ~)
        t            = src.Value;
        sld.Value    = clamp(t, sld.Limits(1), sld.Limits(2));
        val_lbl.Text = sprintf('Value: %.1f', t);
        updateVolume(src, t, viewer, vol);
    end

    function on_range_change(src, ~)
        new_max    = src.Value;
        sld.Limits = [0 new_max];
        spn.Limits = [0 new_max];
        if sld.Value > new_max
            sld.Value = new_max;
            spn.Value = new_max;
            val_lbl.Text = sprintf('Value: %.1f', new_max);
            updateVolume(sld, new_max, viewer, vol);
        end
    end

    function updateVolume(src, t, viewer, vol)
        mask = vol > t;
        mask = bwareaopen(mask, 100);
        data = guidata(src);
        data.mask = mask;
        guidata(src, data);
        viewer.Data = vol .* uint16(mask);
    end

end

% =========================================================================
function v = clamp(v, lo, hi)
    v = max(lo, min(hi, v));
end