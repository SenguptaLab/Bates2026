function annotations = annotate_cells(gfp_volumes)
% ANNOTATE_CELLS  Stage 1 of the pipeline. For each segmented cell volume, the
% user (1) rotates it into a preferred orientation, (2) picks which of the three
% orthogonal mean-projections (XY/XZ/YZ) to keep, and (3) clicks the nucleus.
% The output feeds AVERAGE_CELL_INTENSITY, which only reads .projection and
% .nuc_center (.R and .proj_axis are kept purely as provenance).
%
% External dependency: slanCM (colormap) for the 'thermal' display maps.
%
% Input:
%   gfp_volumes  : 1×N cell array of 3D double, NaN outside cell
%
% Output:
%   annotations  : struct array with fields:
%                    .R            3×3 rotation matrix
%                    .proj_axis    1=XY, 2=XZ, 3=YZ
%                    .nuc_center   [row col] in chosen projection
%                    .projection   2D double mean projection

N           = numel(gfp_volumes);
annotations = struct('R',          cell(1,N), ...
                     'proj_axis',  cell(1,N), ...
                     'nuc_center', cell(1,N), ...
                     'projection', cell(1,N));

for i = 1:N
    fprintf('Cell %d / %d\n', i, N);
    vol = double(gfp_volumes{i});
    ogip = max(vol,[],3);
    imshow(ogip,[])
    [annotations(i).R, vol_rot] = supervise_rotation(vol, i);

    [annotations(i).proj_axis, ...
     annotations(i).nuc_center, ...
     annotations(i).projection] = select_projection(vol_rot, i);
end

end

% =========================================================================
function [R, vol_rot] = supervise_rotation(vol, cell_idx)

    sz          = size(vol);
    mask        = double(~isnan(vol));
    mask_smooth = smooth3(mask, 'gaussian', 5);

    az      = -37.5;
    el      =  30;
    az_step =  5;
    el_step =  5;

    fig = figure('Name', sprintf('Cell %d — Arrow keys=azimuth | Scroll=elevation | Shift+Scroll=azimuth | Enter to confirm', cell_idx), ...
                 'Color', 'k', ...
                 'KeyPressFcn', @on_key, ...
                 'WindowScrollWheelFcn', @on_scroll, ...
                 'Position',[686 446 1159 933]);

    tl    = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax3d  = nexttile(tl, 1);
    ax_xy = nexttile(tl, 2);
    ax_xz = nexttile(tl, 3);
    ax_yz = nexttile(tl, 4);

    for ax = [ax_xy ax_xz ax_yz]
        set(ax, 'Color', 'k'); axis(ax, 'off');
    end

    set(ax3d, 'Color', 'k', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
    axis(ax3d, 'equal'); hold(ax3d, 'on');

    fv = isosurface(mask_smooth, 0.5);
    p  = patch('Parent', ax3d, 'Faces', fv.faces, 'Vertices', fv.vertices, ...
               'FaceColor', [0.3 0.6 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    isonormals(mask_smooth, p);
    set(ax3d, 'Clipping', 'off');
    set(p,    'Clipping', 'off');

    mx = sz(2); my = sz(1); mz = sz(3);

    fill3(ax3d, [1 mx mx 1], [1 1 my my], [mz mz mz mz], ...
          'r', 'FaceAlpha', 0.12, 'EdgeColor', 'r', 'LineWidth', 1.5, 'Clipping', 'off');
    text(ax3d, mx/2, my/2, mz*1.15, 'XY', 'Color', 'r', ...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    quiver3(ax3d, mx/2, my/2, mz*1.3, 0, 0, -mz*0.2, 'r', ...
            'LineWidth', 2, 'MaxHeadSize', 0.5, 'Clipping', 'off');

    fill3(ax3d, [1 mx mx 1], [1 1 1 1], [1 1 mz mz], ...
          'g', 'FaceAlpha', 0.12, 'EdgeColor', 'g', 'LineWidth', 1.5, 'Clipping', 'off');
    text(ax3d, mx/2, -my*0.15, mz/2, 'XZ', 'Color', 'g', ...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    quiver3(ax3d, mx/2, -my*0.3, mz/2, 0, my*0.2, 0, 'g', ...
            'LineWidth', 2, 'MaxHeadSize', 0.5, 'Clipping', 'off');

    fill3(ax3d, [1 1 1 1], [1 my my 1], [1 1 mz mz], ...
          'y', 'FaceAlpha', 0.12, 'EdgeColor', 'y', 'LineWidth', 1.5, 'Clipping', 'off');
    text(ax3d, -mx*0.15, my/2, mz/2, 'YZ', 'Color', 'y', ...
         'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    quiver3(ax3d, -mx*0.3, my/2, mz/2, mx*0.2, 0, 0, 'y', ...
            'LineWidth', 2, 'MaxHeadSize', 0.5, 'Clipping', 'off');

    lighting(ax3d, 'gouraud');
    light('Parent', ax3d, 'Position', [ 1  1  1], 'Style', 'infinite');
    light('Parent', ax3d, 'Position', [-1 -1  0], 'Style', 'infinite');
    axis(ax3d, 'vis3d');
    update_view();

    uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Confirm rotation', ...
              'Units', 'normalized', 'Position', [0.35 0.01 0.3 0.04], ...
              'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'w', ...
              'FontSize', 11, 'Callback', @(~,~) set(fig, 'UserData', 'done'));

    t = timer('ExecutionMode', 'fixedRate', 'Period', 0.3, 'TimerFcn', @update_mips);
    start(t);

    waitfor(fig, 'UserData', 'done');
    stop(t); delete(t);

    view(ax3d, az, el);
    R       = view_to_rotation(ax3d);
    vol_rot = apply_R(vol, R, sz);
    vol_rot = crop_to_content(vol_rot);

    if ishandle(fig); close(fig); end

    % -------------------------------------------------------------------------
    function update_view()
        view(ax3d, az, el);
        title(ax3d, sprintf('Cell %d  |  Az=%.0f  El=%.0f  |  ←→=azimuth  ↑↓=elevation  Scroll=elevation  Shift+Scroll=azimuth', ...
                            cell_idx, az, el), 'Color', 'w', 'FontSize', 8);
        drawnow;
    end

    function on_key(~, evt)
        switch evt.Key
            case 'rightarrow';  az = az + az_step;
            case 'leftarrow';   az = az - az_step;
            case 'uparrow';     el = clamp(el + el_step, -90, 90);
            case 'downarrow';   el = clamp(el - el_step, -90, 90);
            case 'return';      set(fig, 'UserData', 'done'); return;
        end
        update_view();
    end

    function on_scroll(~, evt)
        modifiers = get(fig, 'CurrentModifier');
        if ismember('shift', modifiers)
            % Shift + scroll: azimuth
            az = az - evt.VerticalScrollCount * az_step;
        else
            % Scroll: elevation
            el = clamp(el - evt.VerticalScrollCount * el_step, -90, 90);
        end
        update_view();
    end

    function update_mips(~, ~)
        if ~ishandle(fig); return; end
        try
            view(ax3d, az, el);
            R_cur  = view_to_rotation(ax3d);
            vr     = apply_R(vol, R_cur, sz);
            vr     = crop_to_content(vr);
            draw_mip(ax_xy, mean(vr, 3, 'omitnan'),          'XY', 'r');
            draw_mip(ax_xz, squeeze(mean(vr, 1, 'omitnan')), 'XZ', 'g');
            draw_mip(ax_yz, squeeze(mean(vr, 2, 'omitnan')), 'YZ', 'y');
            drawnow limitrate;
        catch
        end
    end

    function draw_mip(ax, mip, label, col)
        cla(ax);
        imagesc('Parent', ax, 'CData', mip);
        colormap(ax, bone_reverse_steep(256,1));
        axis(ax, 'image', 'off');
        title(ax, label, 'Color', col, 'FontSize', 11, 'FontWeight', 'bold');
        set(ax, 'Color', 'k');
    end

end

% =========================================================================
function [proj_axis, nuc_center, proj] = select_projection(vol_rot, cell_idx)

    mip_xy = mean(vol_rot, 3, 'omitnan');
    mip_xz = squeeze(mean(vol_rot, 1, 'omitnan'));
    mip_yz = squeeze(mean(vol_rot, 2, 'omitnan'));
    mips   = {mip_xy, mip_xz, mip_yz};
    labels = {'XY', 'XZ', 'YZ'};
    colors = {'r',  'g',  'y'};

    fig1 = figure('Name', sprintf('Cell %d — Click projection to keep', cell_idx), ...
                  'Color', 'k','Position',[686 446 1159 933]);

    axes_h = zeros(1, 3);
    for j = 1:3
        axes_h(j) = subplot(1, 3, j, 'Parent', fig1);
        imagesc('Parent', axes_h(j), 'CData', mips{j});
        colormap(axes_h(j), bone_reverse_steep(256,1));
        axis(axes_h(j), 'image', 'off');
        title(axes_h(j), labels{j}, 'Color', colors{j}, 'FontSize', 14, 'FontWeight', 'bold');
    end
    sgtitle(sprintf('Cell %d — Click the projection to keep', cell_idx), 'Color', 'w');

    proj_axis = [];
    set(fig1, 'WindowButtonDownFcn', @on_pick);
    waitfor(fig1, 'UserData', 'picked');
    if ishandle(fig1); close(fig1); end

    proj      = mips{proj_axis};
    proj_sz   = size(proj);
    pos       = round(proj_sz / 2);

    vmin       = min(proj(:), [], 'omitnan');
    vmax       = max(proj(:), [], 'omitnan');
    clim_disp  = [vmin vmax];
    clim_range = vmax - vmin;

    fig2 = figure('Name', sprintf('Cell %d — Click nucleus | Scroll=brightness | Shift+Scroll=contrast | Confirm when done', cell_idx), ...
                  'Color', 'k', 'KeyPressFcn', @on_key2, ...
                  'WindowScrollWheelFcn', @on_scroll2, ...
                  'Position',[686 446 1159 933]);

    ax_proj = subplot(1, 2, 1, 'Parent', fig2);
    ax_info = subplot(1, 2, 2, 'Parent', fig2);

    uicontrol('Parent', fig2, 'Style', 'pushbutton', 'String', 'Confirm nucleus', ...
              'Units', 'normalized', 'Position', [0.35 0.01 0.3 0.05], ...
              'BackgroundColor', [0.2 0.6 0.2], 'ForegroundColor', 'w', ...
              'FontSize', 11, 'Callback', @(~,~) set(fig2, 'UserData', 'done'));

    draw_nuc_view();
    set(fig2, 'WindowButtonDownFcn', @on_click2);
    waitfor(fig2, 'UserData', 'done');
    nuc_center = pos;
    if ishandle(fig2); close(fig2); end

    % -------------------------------------------------------------------------
    function on_pick(~, ~)
        ca = gca;
        for j = 1:3
            if ca == axes_h(j)
                proj_axis = j;
                set(fig1, 'UserData', 'picked');
                return;
            end
        end
    end

    function draw_nuc_view()
        cla(ax_proj);
        imagesc('Parent', ax_proj, 'CData', proj, clim_disp);
        colormap(ax_proj, bone_reverse_steep(256,1));
        axis(ax_proj, 'image', 'off'); hold(ax_proj, 'on');
        plot(ax_proj, pos(2), pos(1), 'r+', 'MarkerSize', 20, 'LineWidth', 2);
        xline(ax_proj, pos(2), 'r');
        yline(ax_proj, pos(1), 'r');
        hold(ax_proj, 'off');
        title(ax_proj, sprintf('%s | Scroll=brightness  Shift+Scroll=contrast', ...
                               labels{proj_axis}), ...
              'Color', colors{proj_axis}, 'FontWeight', 'bold', 'FontSize', 8);

        cla(ax_info); axis(ax_info, 'off');
        text(ax_info, 0.5, 0.75, sprintf('Cell %d', cell_idx), ...
             'Color', 'w', 'FontSize', 14, 'HorizontalAlignment', 'center');
        text(ax_info, 0.5, 0.55, sprintf('Position: [%d, %d]', pos(1), pos(2)), ...
             'Color', 'w', 'FontSize', 11, 'HorizontalAlignment', 'center');
        text(ax_info, 0.5, 0.38, sprintf('CLim: [%.2f  %.2f]', clim_disp(1), clim_disp(2)), ...
             'Color', [0.7 0.7 0.7], 'FontSize', 9, 'HorizontalAlignment', 'center');
        text(ax_info, 0.5, 0.25, 'Click to place nucleus', ...
             'Color', [0.6 0.6 0.6], 'FontSize', 9, 'HorizontalAlignment', 'center');
        text(ax_info, 0.5, 0.12, 'Click Confirm when done', ...
             'Color', [0.5 0.8 0.5], 'FontSize', 9, 'HorizontalAlignment', 'center');
    end

    function on_scroll2(~, evt)
        modifiers = get(fig2, 'CurrentModifier');
        delta     = evt.VerticalScrollCount;
        step      = clim_range * 0.05;
        if ismember('shift', modifiers)
            center    = mean(clim_disp);
            half      = (clim_disp(2) - clim_disp(1)) / 2;
            half      = max(half + delta * step, step);
            clim_disp = [center - half, center + half];
        else
            clim_disp = clim_disp + delta * step;
        end
        draw_nuc_view();
    end

    function on_click2(~, ~)
        if gca ~= ax_proj; return; end
        pt     = get(ax_proj, 'CurrentPoint');
        pos(2) = clamp(round(pt(1,1)), 1, proj_sz(2));
        pos(1) = clamp(round(pt(1,2)), 1, proj_sz(1));
        draw_nuc_view();
    end

    function on_key2(~, evt)
        if strcmp(evt.Key, 'return')
            set(fig2, 'UserData', 'done');
        end
    end

end

% =========================================================================
function R = view_to_rotation(ax)
% Turn the current 3D camera orientation into a rotation matrix, so that the
% orientation the user dialled in on screen becomes the volume's new frame.
% Build an orthonormal camera basis (right/up/forward) via Gram-Schmidt:
    cam_pos    = get(ax, 'CameraPosition');
    cam_target = get(ax, 'CameraTarget');
    cam_up     = get(ax, 'CameraUpVector');

    forward = cam_target - cam_pos;           % where the camera looks
    forward = forward / norm(forward);
    up      = cam_up / norm(cam_up);
    right   = cross(forward, up);             % right = forward x up
    right   = right / norm(right);
    up      = cross(right, forward);          % re-orthogonalise up against right/forward
    up      = up / norm(up);

    % Rows = new axes expressed in old coords. -forward puts the viewing axis
    % last so that the XY projection of the rotated volume matches the screen.
    R = [right; up; -forward];
end

% =========================================================================
function vr = apply_R(v, R, sz0)
% Rotate volume v by R about its centre, into a padded output box large enough
% that no voxel is clipped, with the rotated cell re-centred in that box.
% NB: array subscripts are [row col slice] = [y x z]; imwarp/affine3d work in
% [x y z]. The ([2 1 3]) reorderings swap row<->col to move between the two.
    ctr0    = sz0 / 2;
    ctr_xyz = ctr0([2 1 3]);                   % volume centre in x,y,z order

    % Cell centroid (so we can recentre after rotating).
    mask        = ~isnan(v);
    [x,y,z]     = ind2sub(sz0, find(mask));
    centroid    = mean([x y z], 1);
    cen_xyz     = centroid([2 1 3]);
    cen_rot_xyz = (R * (cen_xyz - ctr_xyz)')' + ctr_xyz;   % where the centroid lands after R

    % Pad by the volume diagonal: guarantees any rotation fits without clipping.
    pad         = ceil(sqrt(sum(sz0.^2)));
    sz_pad      = sz0 + 2*pad;
    pad_ctr_xyz = sz_pad([2 1 3]) / 2;
    trans_xyz   = pad_ctr_xyz - cen_rot_xyz;   % shift that recentres the rotated cell

    % Compose the affine: translate to origin, rotate, translate to padded centre.
    % (affine3d uses row-vector * matrix convention, hence the translations sit in
    %  the bottom row and the product order is T1*Raff*T2.)
    T1   = [1 0 0 0; 0 1 0 0; 0 0 1 0; -ctr_xyz           1];
    Raff = eye(4); Raff(1:3,1:3) = R;
    T2   = [1 0 0 0; 0 1 0 0; 0 0 1 0;  ctr_xyz+trans_xyz  1];
    T    = T1 * Raff * T2;

    ref_in  = imref3d(sz0);
    ref_pad = imref3d(sz_pad);
    vr      = imwarp(v, ref_in, affine3d(T), 'linear', ...
                     'OutputView', ref_pad, 'FillValues', NaN);
end

% =========================================================================
function vr = crop_to_content(vr)
% Trim the NaN padding back off: find the bounding box of non-NaN voxels along
% each axis and keep only that sub-volume.
    mask  = ~isnan(vr);
    any_r = any(any(mask, 2), 3);            % rows that contain any data
    any_c = any(any(mask, 1), 3);            % cols that contain any data
    any_s = squeeze(any(any(mask, 1), 2));   % slices that contain any data

    r = find(any_r); c = find(any_c); s = find(any_s);
    if isempty(r) || isempty(c) || isempty(s); return; end
    vr = vr(r(1):r(end), c(1):c(end), s(1):s(end));
end

% =========================================================================
function v = clamp(v, lo, hi)
    v = max(lo, min(hi, v));
end