function results = average_cell_intensity(annotations)
% AVERAGE_CELL_INTENSITY  Register per-cell 2D projections on a common canvas,
% quantify nuclear/cytoplasmic (N/C) intensity ratios three ways, and build the
% population-average projection + linescan profiles consumed by REPLAY_FIGURE.
%
% Input:
%   annotations : 1xN struct array from ANNOTATE_CELLS, with (at least) fields
%                   .projection  2D mean-projection of one cell (NaN outside cell)
%                   .nuc_center  [row col] of the nucleus in that projection
%
% Output:
%   results : struct consumed by REPLAY_FIGURE and COMBINE_RESULTS_FROM_FILES.
%             Per-cell column vectors (length N): nc_ratios, nc_ratios_mean,
%             col_nc_ratios, nuc_radii. Population fields: avg_proj, avg_crop,
%             avg_smooth, clim, cc_img, canvas_sz, col_profiles, row_profiles,
%             proj_stack, N, annotations, replay.
%
% External dependency: slanCM (colormaps) via the interactive N/C tool.
%
% Coordinate convention that ties this file to REPLAY_FIGURE:
%   Every cell is shifted so its nucleus lands at the canvas CENTRE
%   (row canvas_sz(1)/2, col canvas_sz(2)/2). Everything downstream assumes the
%   nucleus is at the canvas centre, so the final crop is kept symmetric about
%   that centre in BOTH axes.

N = numel(annotations);

% --- Canvas sizing -------------------------------------------------------
% Canvas is 2x the largest projection so that, after we recentre each cell's
% nucleus to the middle, no cell can be pushed off the edge.
proj_sizes = cellfun(@(p) size(p), {annotations.projection}, 'UniformOutput', false);
max_sz     = max(cat(1, proj_sizes{:}), [], 1);
canvas_sz  = 2 * max_sz;

% --- Preallocation -------------------------------------------------------
proj_stack     = NaN([canvas_sz N]);   % every cell's registered projection, stacked in 3rd dim
nc_ratios      = NaN(N, 1);            % N/C ratio, radial method (nucleus disc vs brightest cyto ring)
nc_ratios_adj  = NaN(N, 1);            % N/C ratio, adjacent method (nucleus vs brightest CONNECTED cyto patch)
nc_ratios_mean = NaN(N, 1);            % N/C ratio, simple method (mean-in vs mean-out)
nuc_radii      = NaN(N, 1);            % nucleus radius (px) set interactively
col_nc_ratios  = NaN(N, 1);           % N/C ratio, 1D horizontal-linescan method
col_profiles   = NaN(canvas_sz(2), N); % per-cell column-mean profile (intensity vs column)
row_profiles   = NaN(canvas_sz(2), N); % per-cell single-row linescan through the nucleus
ls_left_cols   = NaN(N, 2);            % per-cell LEFT cytoplasm flank window [start end], stored
ls_right_cols  = NaN(N, 2);            % as OFFSETS from the nucleus column (registration-invariant)
nuc_sig_frac   = NaN(N, 1);            % share of total (min-subtracted) cell signal inside the nucleus [0..1]
gini_coef      = NaN(N, 1);            % Gini coefficient of cell pixel intensities (concentration/inequality)
entropy_coef   = NaN(N, 1);            % 1 - normalized Shannon entropy of cell pixel intensities (concentration)
ls_gini_coef   = NaN(N, 1);            % Gini coefficient restricted to the linescan row only
ls_entropy_coef = NaN(N, 1);           % 1 - normalized Shannon entropy restricted to the linescan row only

% Canvas-centre coordinates. Nucleus is registered to (rc_mid, cc_mid).
rc_mid = canvas_sz(1) / 2;
cc_mid = canvas_sz(2) / 2;

% --- Radial-profile setup ------------------------------------------------
% Radial profile = mean intensity in concentric 1-px rings about the nucleus
% centre (which is the canvas centre after registration). Because the centre is
% fixed, the pixel->ring assignment is identical for every cell, so precompute it
% once. Ring k (k = 1..r_max) covers distances [k-1, k) px from the centre. Radius
% is a distance-from-nucleus index -> registration-invariant, so it pools cleanly
% across differently sized canvases in COMBINE_RESULTS_FROM_FILES (pad the shorter
% profiles at their FAR end, no re-centring needed).
[CC_grid, RR_grid] = meshgrid(1:canvas_sz(2), 1:canvas_sz(1));
dist_grid = sqrt((RR_grid - rc_mid).^2 + (CC_grid - cc_mid).^2);
r_max     = floor(min(canvas_sz) / 2);          % keep rings within the canvas on all sides
ring_bin  = min(r_max, floor(dist_grid) + 1);   % ring index for every pixel, clamped to r_max
radial_profiles = NaN(r_max, N);                % per-cell mean intensity vs radius (bin centres below)
radial_sums     = NaN(r_max, N);                % per-cell total intensity vs radius (empty rings -> 0)
radial_r_px     = (0:r_max-1)' + 0.5;           % radius (px) at each ring centre
radial_half_px  = NaN(N, 1);                    % radius enclosing 50% of the cell's signal (median radius)
radial_halfmax_px = NaN(N, 1);                  % radius where the profile first reaches half its own max

for i = 1:N
    proj = double(annotations(i).projection);
    mask = ~isnan(proj);

    % --- 1. Min-max normalize within cell ---
    % Done first because the scaling is per-pixel and therefore commutes with the
    % translation in step 2; order doesn't affect the result.
    mu  = min(proj(mask));
    rng = max(proj(mask)) - mu;
    if rng > 0
        proj = (proj - mu) / rng;
    end
    proj(~mask) = NaN;

    % --- 2. Register: place the cell on the canvas with its nucleus at the centre ---
    % place_on_canvas_2d lands the given point (here, the nucleus) directly on the
    % canvas centre, so no pre-shift is needed. This is a single clip-safe
    % translation onto the 2x-larger canvas. It replaces an earlier circshift-based
    % recentering: circshift wraps pixels around a same-size array and the wrapped
    % band then had to be blanked, which SILENTLY DISCARDED real cell pixels on the
    % side away from the nucleus (up to ~half the cell when the nucleus sat near a
    % projection edge). Direct placement can't wrap, so every cell pixel is kept.
    canvas            = place_on_canvas_2d(proj, canvas_sz, annotations(i).nuc_center);
    proj_stack(:,:,i) = canvas;

    % --- 3. Column mean profile ---
    col_profiles(:,i) = mean(canvas, 1, 'omitnan');

    % --- 4. Horizontal linescan through nucleus row ---
    nuc_row           = round(canvas_sz(1) / 2);
    row_profiles(:,i) = canvas(nuc_row, :)';

    % --- 4b. Radial profiles about the nucleus centre ---
    % Mean intensity in each 1-px ring (empty rings -> NaN), and total intensity
    % (sum) in each ring (empty rings -> 0). The sum grows with ring circumference
    % even for uniform signal; its cumulative (built in REPLAY_FIGURE) gives the
    % fraction of the cell's signal within each radius -- a spatial CDF whose SHAPE
    % cleanly separates nuclear-included (rises fast near centre) from excluded
    % (rises late), which the whole-cell concentration scalars barely register.
    vals_i  = canvas(:);
    ok_i    = ~isnan(vals_i);
    radial_profiles(:, i) = accumarray(ring_bin(ok_i), vals_i(ok_i), [r_max 1], @mean, NaN);
    radial_sums(:, i)     = accumarray(ring_bin(ok_i), vals_i(ok_i), [r_max 1], @sum, 0);
    % Median radius: interpolated radius enclosing 50% of the cell's total signal.
    % A single scalar read-out of the cumulative radial CDF -- small = signal
    % concentrated near the nucleus centre, large = pushed to the periphery -- so
    % it discriminates nuclear inclusion from exclusion symmetrically.
    radial_half_px(i)     = radial_frac_radius(radial_sums(:, i), radial_r_px, 0.5);
    % Half-max radius: the radius at which this cell's OWN radial profile first
    % reaches half of its own maximum, scanning outward from the centre. Unlike
    % r50 (which is cumulative, a property of the whole signal distribution), this
    % is a shape feature of the profile curve itself -- a width/decay scale. It is
    % direction-agnostic: for a nuclear-INCLUDED profile (peak at the centre) the
    % first crossing is the DECAY radius on the way down; for a nuclear-EXCLUDED
    % profile (dip at the centre, rising outward) it is the RISE radius on the way
    % up. Level is half the profile max; see the helper for the baseline caveat.
    radial_halfmax_px(i)  = radial_half_max_radius(radial_profiles(:, i), radial_r_px);

    % --- 5. Interactive N/C ratio annotation (radial + mean) ---
    ctr_r = canvas_sz(1) / 2;
    ctr_c = canvas_sz(2) / 2;
    [nc_ratios(i), nc_ratios_mean(i), nuc_radii(i), nuc_sig_frac(i), gini_coef(i), entropy_coef(i), ...
        nc_ratios_adj(i)] = annotate_nc_ratio(canvas, ctr_r, ctr_c, i);

    % --- 6. Linescan-based N/C ratio ---
    % Work along the single horizontal row through the nucleus. Nucleus signal =
    % mean over the nuclear disc width; cytoplasm signal = the pixels in BOTH
    % flanking peaks (the membrane/cytoplasm shows up as bright shoulders on each
    % side of the nucleus).
    nuc_col      = round(ctr_c);              % nucleus column = canvas centre
    linescan_row = canvas(nuc_row, :);        % the 1D trace through the nucleus
    nuc_r_px     = round(nuc_radii(i));       % nucleus half-width from the interactive step

    % Nuclear window: [-r, +r] about the centre.
    nuc_cols    = max(1, nuc_col - nuc_r_px) : min(canvas_sz(2), nuc_col + nuc_r_px);
    nuc_ls_vals = linescan_row(nuc_cols);
    nuc_ls_mean = mean(nuc_ls_vals, 'omitnan');

    % LEFT flank: search a band 1r..3r left of the nucleus, find its brightest
    % column (the cytoplasm/membrane peak), then average a +/- r window about it,
    % clipped so it never overlaps the nuclear window.
    left_search_start = max(1, nuc_col - 3*nuc_r_px);
    left_search_end   = max(1, nuc_col - nuc_r_px - 1);
    left_search_vals  = linescan_row(left_search_start:left_search_end);
    [~, left_max_idx] = max(left_search_vals);
    left_max_col      = left_search_start + left_max_idx - 1;

    left_flank_start = max(1, left_max_col - nuc_r_px);
    left_flank_end   = min(canvas_sz(2), left_max_col + nuc_r_px);
    left_flank_end   = min(left_flank_end, nuc_col - nuc_r_px - 1);   % keep off the nucleus
    left_vals        = linescan_row(left_flank_start:left_flank_end);

    % RIGHT flank: mirror image of the left-flank logic.
    right_search_start = min(canvas_sz(2), nuc_col + nuc_r_px + 1);
    right_search_end   = min(canvas_sz(2), nuc_col + 3*nuc_r_px);
    right_search_vals  = linescan_row(right_search_start:right_search_end);
    [~, right_max_idx] = max(right_search_vals);
    right_max_col      = right_search_start + right_max_idx - 1;

    right_flank_start = max(1, right_max_col - nuc_r_px);
    right_flank_end   = min(canvas_sz(2), right_max_col + nuc_r_px);
    right_flank_start = max(right_flank_start, nuc_col + nuc_r_px + 1); % keep off the nucleus
    right_vals        = linescan_row(right_flank_start:right_flank_end);

    % Cytoplasm reference = pixels from BOTH flank windows pooled together
    % (previously only the brighter side was used). Pooling the raw values lets a
    % wider/better-sampled flank contribute proportionally; if one side's window is
    % empty (nucleus near an edge) it simply drops out and the other side is used.
    cyto_ls_vals     = [left_vals, right_vals];
    cyto_ls_mean     = mean(cyto_ls_vals, 'omitnan');
    col_nc_ratios(i) = nuc_ls_mean / cyto_ls_mean;

    % Record the columns each flank used, as OFFSETS from the nucleus column. The
    % nucleus is always registered to the canvas centre, so these offsets are
    % independent of canvas size and survive pooling in COMBINE_RESULTS_FROM_FILES
    % unchanged. To recover absolute columns on any canvas (e.g. for shading the
    % individual traces in REPLAY_FIGURE): nucleus_col + offset, where nucleus_col
    % is that struct's canvas_sz(2)/2. Empty windows stay NaN.
    if left_flank_start <= left_flank_end
        ls_left_cols(i, :)  = [left_flank_start,  left_flank_end]  - nuc_col;
    end
    if right_flank_start <= right_flank_end
        ls_right_cols(i, :) = [right_flank_start, right_flank_end] - nuc_col;
    end

    % Concentration metrics restricted to the LINESCAN ROW only (not the whole
    % cell). The whole-cell versions (gini_coef/entropy_coef below) are diluted by
    % everything off this row -- background noise, puncta, cell-shape irregularity
    % elsewhere -- which acts as a roughly constant baseline "concentration"
    % present in every cell regardless of nuclear localization, so a real nuclear
    % effect has to move the needle on top of that baseline. Restricting to the
    % single row that actually crosses the nucleus and both flanks removes most of
    % that off-axis dilution, at the cost of a much smaller per-cell sample (tens
    % to ~100 pixels vs. thousands), so individual-cell values are noisier; average
    % across cells in a group as usual.
    ls_vals            = linescan_row(~isnan(linescan_row));
    ls_gini_coef(i)     = gini(ls_vals);
    ls_entropy_coef(i)  = shannon_conc(ls_vals);

    fprintf(['Cell %d: N/C (radial) = %.3f  |  N/C (adj) = %.3f  |  N/C (mean) = %.3f  |  N/C (linescan) = %.3f  ' ...
             '|  nuc frac = %.3f  |  Gini = %.3f  |  entropy = %.3f  |  r50 = %.1f px  |  r_halfmax = %.1f px  (nuc_r = %d px)\n'], ...
            i, nc_ratios(i), nc_ratios_adj(i), nc_ratios_mean(i), col_nc_ratios(i), ...
            nuc_sig_frac(i), gini_coef(i), entropy_coef(i), radial_half_px(i), ...
            radial_halfmax_px(i), round(nuc_radii(i)));
end

% --- 7. Average only well-covered pixels ---
% Average across cells, but null out any pixel that fewer than min_coverage*N
% cells contributed to, so thin edge regions don't create noisy averages.
min_coverage = 0.5;
n_per_pixel  = sum(~isnan(proj_stack), 3);
avg_proj     = mean(proj_stack, 3, 'omitnan');
avg_proj(n_per_pixel < min_coverage * N) = NaN;

% --- 8. Crop symmetrically about the nucleus centre (both axes) ---
% The nucleus sits at (rc_mid, cc_mid). We crop the same distance either side of
% that centre in rows AND columns so the nucleus stays dead-centre in the output
% image. REPLAY_FIGURE relies on this: it draws the linescan row at img_rows/2 and
% the nucleus column at cc_img, both of which are only correct for a centred crop.
% (Previously only the columns were centred, so the drawn linescan row could miss
%  the nucleus; matching COMBINE_RESULTS_FROM_FILES here fixes that.)
mask_avg = ~isnan(avg_proj);
r        = find(any(mask_avg, 2));   % occupied rows
c        = find(any(mask_avg, 1));   % occupied columns
margin   = 20;

half_h = max(rc_mid - r(1), r(end) - rc_mid);   % furthest occupied row from centre
half_w = max(cc_mid - c(1), c(end) - cc_mid);   % furthest occupied col from centre

r1 = max(1,                round(rc_mid - half_h) - margin);
r2 = min(size(avg_proj,1), round(rc_mid + half_h) + margin);
c1 = max(1,                round(cc_mid - half_w) - margin);
c2 = min(size(avg_proj,2), round(cc_mid + half_w) + margin);

avg_crop  = avg_proj(r1:r2, c1:c2);
mask_crop = ~isnan(avg_crop);
cc_img    = cc_mid - c1 + 1;          % nucleus column expressed in cropped coords
rr_img    = rc_mid - r1 + 1;          % nucleus row expressed in cropped coords (mirror of cc_img)
                                      % Stored explicitly rather than assumed to be img_rows/2:
                                      % if the crop clips against the canvas edge the crop is no
                                      % longer symmetric about the centre, so REPLAY_FIGURE must be
                                      % told the true nucleus row instead of guessing the midpoint.

% --- 9. Fill NaN boundary with zeros, then lightly smooth for display ---
avg_filled             = avg_crop;
avg_filled(~mask_crop) = 0;
avg_smooth             = imgaussfilt(avg_filled, 1);

% Display limits. Lower bound is fixed at 0 because REPLAY_FIGURE forces
% clim_display = [0, clim(2)]; we only need the 99th-percentile upper bound here.
vals = avg_smooth(mask_crop);
clim = [0, prctile(vals, 99)];

% --- 10. Build results struct ---
% Field set below is the shared "contract": REPLAY_FIGURE reads avg_smooth, clim,
% the four per-cell vectors, col/row_profiles, cc_img, canvas_sz, N;
% COMBINE_RESULTS_FROM_FILES additionally reads avg_proj, proj_stack, annotations.
results.avg_proj       = avg_proj;      % full averaged canvas (used by combine)
results.avg_crop       = avg_crop;      % cropped, unsmoothed average (kept as raw reference)
results.avg_smooth     = avg_smooth;    % cropped, smoothed average (displayed by replay)
results.clim           = clim;
results.nc_ratios      = nc_ratios;
results.nc_ratios_adj  = nc_ratios_adj;   % adjacent-cytoplasm radial method
results.nc_ratios_mean = nc_ratios_mean;
results.col_nc_ratios  = col_nc_ratios;
results.nuc_sig_frac   = nuc_sig_frac;   % share of total cell signal inside the nucleus [0..1]
results.gini_coef      = gini_coef;      % Gini coefficient of cell pixel intensities
results.entropy_coef   = entropy_coef;   % 1 - normalized Shannon entropy of cell pixel intensities
results.ls_gini_coef    = ls_gini_coef;     % Gini restricted to the linescan row
results.ls_entropy_coef = ls_entropy_coef;  % 1 - normalized Shannon entropy restricted to the linescan row
results.nuc_radii      = nuc_radii;
results.col_profiles   = col_profiles;
results.row_profiles   = row_profiles;
results.radial_profiles = radial_profiles;  % per-cell mean intensity vs radius from nucleus centre
results.radial_sums     = radial_sums;       % per-cell total intensity vs radius (for the cumulative radial CDF)
results.radial_half_px  = radial_half_px;    % radius enclosing 50% of signal (median radius; CDF x@0.5)
results.radial_halfmax_px = radial_halfmax_px; % radius where the profile first reaches half its own max
results.radial_r_px     = radial_r_px;      % radius (px) at each ring centre (row-matches radial_profiles)
results.ls_left_cols   = ls_left_cols;   % per-cell LEFT flank window [start end], offset from nucleus col
results.ls_right_cols  = ls_right_cols;  % per-cell RIGHT flank window [start end], offset from nucleus col
results.cc_img         = cc_img;
results.rr_img         = rr_img;
results.canvas_sz      = canvas_sz;
results.proj_stack     = proj_stack;
results.N              = N;
results.annotations    = annotations;

% Convenience handle to regenerate every figure from the saved struct. The
% closure captures `results` as it stands here (without .replay itself), which is
% all REPLAY_FIGURE needs. NOTE: previously mis-typed as `cresults.replay`, so the
% handle was silently dropped -- now correctly attached to `results`.
results.replay = @() replay_figure(results);

% --- 11. Generate plots now ---
replay_figure(results);

end

% =========================================================================
function [nc_ratio, nc_ratio_mean, nuc_r_out, nuc_sig_frac, gini_coef, entropy_coef, nc_ratio_adj] = annotate_nc_ratio(canvas, ctr_r, ctr_c, cell_idx)
% Interactive: user scrolls to size a circle over the nucleus, presses Enter.
% Returns two N/C ratios, the chosen radius, and three concentration metrics:
%   nc_ratio      (radial) = mean(nucleus disc) / mean(brightest cyto pixels,
%                            matched in pixel count to the nucleus)
%   nc_ratio_mean (simple) = mean(nucleus disc) / mean(all cytoplasm pixels)
%   nuc_sig_frac  = sum(nucleus) / sum(whole cell) -- share of the cell's total
%                   (min-subtracted) signal inside the nucleus. A ratio of SUMS
%                   rather than means, so it moves strongly as signal redistributes
%                   into the nucleus and is bounded [0,1].
%   gini_coef     = Gini coefficient of ALL cell-pixel intensities (0 = uniform,
%                   ->1 = concentrated). Mask-independent, so it does not depend on
%                   how the nucleus circle was drawn.
%   entropy_coef  = 1 - (normalized Shannon entropy) of ALL cell-pixel intensities,
%                   treated as a probability distribution over pixel LOCATIONS
%                   (p_i = pixel_i / total signal). 0 = signal spread perfectly
%                   uniformly across every pixel (max entropy), ->1 = all signal
%                   piled into a single pixel (min entropy). Same [0,1] direction
%                   as gini_coef (higher = more concentrated) for easy comparison,
%                   but a mathematically distinct measure (Shannon information
%                   rather than the Lorenz-curve-based Gini) and, like Gini,
%                   independent of the nucleus circle.

    fig   = figure('Name', sprintf('Cell %d — Scroll to resize, Enter to confirm', cell_idx), ...
                   'Color', 'k', 'KeyPressFcn', @on_key, ...
                   'WindowScrollWheelFcn', @on_scroll);
    ax    = axes('Parent', fig, 'Color', 'k');
    nuc_r = 10;                                  % initial nucleus radius (px)

    draw_circle();
    waitfor(fig, 'UserData', 'done');            % block until user confirms

    % Split the cell into nucleus (inside circle) and cytoplasm (outside) by
    % radial distance from the marked centre; ignore NaN (outside-cell) pixels.
    [cc, rr]  = meshgrid(1:size(canvas,2), 1:size(canvas,1));
    dist_px   = sqrt((rr - ctr_r).^2 + (cc - ctr_c).^2);
    mask_c    = ~isnan(canvas);
    nuc_mask  = dist_px <= nuc_r & mask_c;
    cyto_mask = dist_px >  nuc_r & mask_c;

    if any(nuc_mask(:)) && any(cyto_mask(:))
        nuc_mean      = mean(canvas(nuc_mask));
        % Radial method: compare the nucleus to the brightest cytoplasm pixels,
        % taking exactly as many cyto pixels as there are nucleus pixels. This
        % targets the bright peri-nuclear cytoplasm rather than diluting with the
        % dim cell periphery.
        N_px          = sum(nuc_mask(:));
        cyto_vals     = sort(canvas(cyto_mask), 'descend');
        N_px          = min(N_px, numel(cyto_vals));
        cyto_mean_top = mean(cyto_vals(1:N_px));
        nc_ratio      = nuc_mean / cyto_mean_top;
        % Simple method: nucleus vs. ALL cytoplasm.
        cyto_mean_all = mean(canvas(cyto_mask));
        nc_ratio_mean = nuc_mean / cyto_mean_all;
        % Adjacent-cytoplasm radial method: like the radial method, but the
        % cytoplasm comparison must be a single CONNECTED (adjacent) group of
        % N_px pixels rather than the brightest pixels scattered anywhere. This
        % targets a coherent bright cytoplasmic patch (a real local accumulation)
        % instead of scattered bright noise, which the plain top-N method can pick
        % up. Found by greedy region-growing from the brightest cytoplasm seeds.
        adj_cyto_mean = brightest_adjacent_cyto_mean(canvas, cyto_mask, nuc_r);
        nc_ratio_adj  = nuc_mean / adj_cyto_mean;
    else
        nc_ratio      = NaN;
        nc_ratio_mean = NaN;
        nc_ratio_adj  = NaN;
    end
    nuc_r_out = nuc_r;

    % --- Concentration metrics ---
    % Nuclear signal fraction: share of the cell's total signal that sits inside
    % the nucleus disc. Because canvas is min-max normalised per cell, its floor is
    % 0 (an implicit background subtraction) and the ratio of sums is scale-free.
    % Tracks signal REDISTRIBUTING into the nucleus far more strongly than a ratio
    % of means (a mean-out is diluted by the large, mostly-unchanged periphery).
    cell_vals  = canvas(mask_c);
    tot_signal = sum(cell_vals);
    if tot_signal > 0
        nuc_sig_frac = sum(canvas(nuc_mask)) / tot_signal;
    else
        nuc_sig_frac = NaN;
    end

    % Gini coefficient over every cell pixel: how unequally signal is distributed
    % across the cell, independent of the nucleus circle. A sharp concentration
    % (into the nucleus or anywhere) drives it up.
    gini_coef = gini(cell_vals);

    % Shannon entropy over every cell pixel, treated as a spatial probability
    % distribution (not an intensity histogram): a signal spread evenly across
    % every pixel has maximum entropy; a signal piled into few pixels has low
    % entropy. Normalized by log2(n_pixels) so cells of different sizes are
    % comparable, then flipped (1 - .) so higher = more concentrated, matching
    % gini_coef's direction.
    entropy_coef = shannon_conc(cell_vals);

    if ishandle(fig); close(fig); end

    function draw_circle()
        cla(ax);
        vmin = min(canvas(:), [], 'omitnan');
        vmax = max(canvas(:), [], 'omitnan');
        cd   = (canvas - vmin) / (vmax - vmin);
        cd(isnan(cd)) = 0;
        imagesc('Parent', ax, 'CData', cd);
        colormap(ax, slanCM('thermal')); axis(ax, 'image', 'off');
        hold(ax, 'on');
        th = linspace(0, 2*pi, 200);
        plot(ax, ctr_c + nuc_r*cos(th), ctr_r + nuc_r*sin(th), 'w-', 'LineWidth', 2);
        plot(ax, ctr_c, ctr_r, 'w+', 'MarkerSize', 12, 'LineWidth', 2);
        hold(ax, 'off');
        title(ax, sprintf('Cell %d  |  Scroll to resize (r=%d px)  |  Enter to confirm', ...
                          cell_idx, round(nuc_r)), 'Color', 'w');
    end

    function on_scroll(~, evt)
        nuc_r = max(1, nuc_r - evt.VerticalScrollCount * 2);
        draw_circle();
    end

    function on_key(~, evt)
        if strcmp(evt.Key, 'return')
            set(fig, 'UserData', 'done');
        end
    end

end

% =========================================================================
function canvas = place_on_canvas_2d(proj, canvas_sz, nuc)
% Paste `proj` onto a blank canvas so that the point `nuc` (in proj coords)
% lands on the canvas centre. Any part of proj that would fall outside the
% canvas is clipped; both the destination and source index ranges are adjusted
% together so the copied blocks always match in size.
    canvas     = NaN(canvas_sz);
    proj_sz    = size(proj);
    canvas_ctr = round(canvas_sz / 2);
    offset     = canvas_ctr - nuc;                 % translation proj -> canvas

    % Desired destination block (may spill past the canvas edges).
    r1 = 1 + offset(1);  r2 = proj_sz(1) + offset(1);
    c1 = 1 + offset(2);  c2 = proj_sz(2) + offset(2);

    % If the block starts before pixel 1, skip that many source pixels.
    src_r1 = 1 + max(0, 1 - r1);
    src_c1 = 1 + max(0, 1 - c1);

    % Clip destination to the canvas bounds.
    r1 = max(1, r1);  r2 = min(canvas_sz(1), r2);
    c1 = max(1, c1);  c2 = min(canvas_sz(2), c2);

    % Matching source end indices (same width/height as the clipped destination).
    src_r2 = src_r1 + (r2 - r1);
    src_c2 = src_c1 + (c2 - c1);

    canvas(r1:r2, c1:c2) = proj(src_r1:src_r2, src_c1:src_c2);
end

% =========================================================================
function m = brightest_adjacent_cyto_mean(canvas, cyto_mask, nuc_r)
% Mean intensity of the brightest CONNECTED (adjacent) group of cytoplasmic
% pixels the same size as the nucleus. Realised as the brightest fully-
% cytoplasmic DISK of radius nuc_r (same disk as the nucleus, so the group has
% the same pixel count): slide that disk over the image and take the position
% with the highest mean. A disk is a compact adjacent patch, so -- unlike the
% plain "brightest N pixels anywhere" radial method -- this rejects scattered
% bright noise pixels and only scores a genuine coherent bright cytoplasmic
% accumulation. Computed by convolution (fast, deterministic).
    if nuc_r < 1 || ~any(cyto_mask(:))
        m = NaN;
        return;
    end
    rd = max(1, round(nuc_r));
    [kx, ky] = meshgrid(-rd:rd, -rd:rd);
    D    = double((kx.^2 + ky.^2) <= rd^2);      % disk, same radius as the nucleus
    area = sum(D(:));

    cf              = canvas;
    cf(~cyto_mask)  = 0;                          % zero out non-cytoplasm before summing
    cyto_sum        = conv2(cf, D, 'same');       % sum of cyto intensity under the disk
    cyto_cnt        = conv2(double(cyto_mask), D, 'same');  % # cyto pixels under the disk
    local_mean      = cyto_sum ./ max(cyto_cnt, 1);

    full = cyto_cnt >= area - 0.5;                % disk entirely within cytoplasm
    if any(full(:))
        m = max(local_mean(full));
    else
        % Thin cytoplasm: no full-size disk fits. Fall back to the most-covered
        % disk position(s) so the metric is still defined (approximate size).
        cand = cyto_cnt >= max(cyto_cnt(:)) - 0.5;
        m = max(local_mean(cand));
    end
end

% =========================================================================
function g = gini(x)
% GINI  Gini coefficient of a set of values (measure of concentration/inequality).
% 0 = perfectly uniform, approaches 1 as all the signal piles into a few pixels.
% NaNs are dropped and values are floored at 0 (Gini requires non-negative input).
    x = x(:);
    x = x(~isnan(x));
    x = max(0, x);
    n = numel(x);
    s = sum(x);
    if n == 0 || s <= 0
        g = NaN;
        return;
    end
    x   = sort(x);
    idx = (1:n)';
    % Sorted-order formula: G = (2*Σ i*x_i)/(n*Σx) - (n+1)/n.
    g = (2 * sum(idx .* x)) / (n * s) - (n + 1) / n;
end

% =========================================================================
function h = shannon_conc(x)
% SHANNON_CONC  1 - normalized Shannon entropy of a set of values, treated as a
% probability distribution over pixel LOCATIONS (p_i = x_i / sum(x)), NOT an
% intensity histogram. 0 = signal spread perfectly evenly across every pixel
% (maximum entropy); ->1 = all signal concentrated into a single pixel (minimum
% entropy). Same direction and [0,1] range as GINI, so the two are directly
% comparable, but Shannon entropy and the Gini coefficient are mathematically
% distinct summaries of "how concentrated" a distribution is, so agreement
% between them is a stronger signal than either alone.
% NaNs are dropped and values are floored at 0 (entropy requires non-negative,
% and p_i = 0 contributes 0 to the sum by convention, i.e. 0*log(0) := 0).
    x = x(:);
    x = x(~isnan(x));
    x = max(0, x);
    n = numel(x);
    s = sum(x);
    if n <= 1 || s <= 0
        h = NaN;
        return;
    end
    p        = x / s;
    p_nz     = p(p > 0);                    % drop zeros: 0*log2(0) := 0
    H        = -sum(p_nz .* log2(p_nz));    % Shannon entropy, in bits
    H_max    = log2(n);                     % entropy of a perfectly uniform distribution
    h        = 1 - H / H_max;
end

% =========================================================================
function r = radial_frac_radius(ring_sum, r_px, level)
% Interpolated radius at which the cumulative radial signal first reaches `level`
% (e.g. 0.5 -> radius enclosing half the total signal, i.e. the median radius).
% ring_sum = per-ring total intensity (radial_sums column); r_px = ring-centre
% radii. Empty rings count as 0 signal. Returns NaN if there is no signal.
    ring_sum = ring_sum(:);
    ring_sum(isnan(ring_sum)) = 0;
    tot = sum(ring_sum);
    if tot <= 0
        r = NaN;
        return;
    end
    cf = cumsum(ring_sum) / tot;
    k  = find(cf >= level, 1, 'first');
    if isempty(k)
        r = r_px(end);
    elseif k == 1
        r = r_px(1);
    else
        c0 = cf(k-1);  c1 = cf(k);
        if c1 > c0
            r = r_px(k-1) + (level - c0) / (c1 - c0) * (r_px(k) - r_px(k-1));
        else
            r = r_px(k);
        end
    end
end