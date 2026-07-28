function replay_figure(results, show_extra_figs)
% REPLAY_FIGURE  Regenerate every figure from a results struct produced by
% AVERAGE_CELL_INTENSITY or COMBINE_RESULTS_FROM_FILES. Pure display: it reads
% the struct and draws, computing all crop/flank geometry on the fly (so nothing
% here can go stale relative to the stored struct).
%
% Reads exactly these fields:
%   avg_smooth, clim, nc_ratios, nc_ratios_mean, col_nc_ratios, nuc_radii,
%   col_profiles, row_profiles, cc_img, canvas_sz, N
% Note: only clim(2) is used — the displayed lower limit is forced to 0.
% Also reads avg_crop if present (used to round the mask edge in Figure 1); older
% structs without it fall back to an approximate mask derived from avg_smooth.
%
% show_extra_figs (optional, default true): when false, ONLY Figure 1 (the average
% projection + linescan panel) is drawn and Figures 2-7 are skipped. REPLAY_FIGURE
% always draws Figure 1.
%
% Figures:
%   1  average projection + segment-coloured horizontal linescan (mean±SEM)
%   2  all three N/C methods as jittered points with group means   [extra]
%   3-5 N/C distributions (radial / simple mean / linescan)        [extra]
%   6  per-cell column profiles + mean                              [extra]
%   7  per-cell horizontal linescans + mean±SEM                     [extra]
%   8  concentration metrics (nuc frac / Gini / entropy / ls variants) [extra]
%   9  radial intensity profile about the nucleus centre + mean±SEM [extra]
%   10 radial signal SUM per ring + mean                            [extra]
%   11 cumulative radial signal distribution (spatial CDF) + mean   [extra]
%

if nargin < 2 || isempty(show_extra_figs); show_extra_figs = true; end

avg_smooth     = results.avg_smooth;
clim           = results.clim;
nc_ratios      = results.nc_ratios;
nc_ratios_mean = results.nc_ratios_mean;
col_nc_ratios  = results.col_nc_ratios;
nuc_radii      = results.nuc_radii;
col_profiles   = results.col_profiles;
row_profiles   = results.row_profiles;
cc_img         = results.cc_img;
canvas_sz      = results.canvas_sz;
% New concentration metrics + radial profile (optional: older structs lack them,
% in which case the corresponding extra figures are skipped).
getf = @(name) getfield_or_empty(results, name);
nuc_sig_frac    = getf('nuc_sig_frac');
gini_coef       = getf('gini_coef');
entropy_coef    = getf('entropy_coef');
ls_gini_coef    = getf('ls_gini_coef');
ls_entropy_coef = getf('ls_entropy_coef');
radial_profiles = getf('radial_profiles');
radial_sums     = getf('radial_sums');
radial_half_px  = getf('radial_half_px');
radial_halfmax_px = getf('radial_halfmax_px');
radial_r_px     = getf('radial_r_px');
% Nucleus row within the crop. Newer structs store it explicitly (mirror of
% cc_img); for older structs fall back to the image midpoint, which is exact
% whenever the row-crop wasn't clipped against the canvas edge.
if isfield(results, 'rr_img')
    rr_img = results.rr_img;
else
    rr_img = round(size(avg_smooth, 1) / 2);
end
N              = results.N;
img_rows       = size(avg_smooth, 1);
img_cols       = size(avg_smooth, 2);
cc_mid         = canvas_sz(2) / 2;
cmap           = lines(N);
x_full         = reshape(1:canvas_sz(2), 1, []);
x_crop_img     = reshape(1:img_cols, 1, []);

col_mean_full = mean(col_profiles, 2, 'omitnan')';
col_sem_full  = (std(col_profiles, 0, 2, 'omitnan') / sqrt(N))';
row_mean_full = mean(row_profiles, 2, 'omitnan')';
row_sem_full  = (std(row_profiles, 0, 2, 'omitnan') / sqrt(N))';
valid_row     = ~isnan(row_mean_full);
valid_col     = ~isnan(col_mean_full);

nuc_row_img = rr_img;   % true nucleus row in cropped coords (see rr_img above)

% --- Crop row profiles to image coords ---
% row_profiles are stored on the FULL canvas (one row per canvas column). The
% displayed image is a crop of that canvas. cc_img is the nucleus column within
% the crop and cc_mid is the nucleus column on the full canvas, so their
% difference recovers where the crop begins on the full canvas. We then slice the
% profiles to the same window the image occupies, keeping the two panels aligned.
c_start = round(cc_mid - cc_img + 1);
c_end   = c_start + img_cols - 1;
c_end   = min(c_end, size(row_profiles, 1));
c_start = max(1, c_start);

row_prof_crop  = row_profiles(c_start:c_end, :)';
row_mean_img   = mean(row_prof_crop, 1, 'omitnan');
row_sem_img    = std(row_prof_crop, 0, 1, 'omitnan') / sqrt(N);
valid_mean_img = ~isnan(row_mean_img);

% --- Nucleus boundaries in image coords ---
mean_nuc_r_px = mean(nuc_radii, 'omitnan');
nuc_col       = round(cc_img);
nuc_r_px      = round(mean_nuc_r_px);
nuc_left      = cc_img - mean_nuc_r_px;
nuc_right     = cc_img + mean_nuc_r_px;

% --- Recompute flank positions from the POPULATION-MEAN trace ---
% Same left/right flank-finding logic as the per-cell linescan in
% AVERAGE_CELL_INTENSITY, but applied to the mean image trace, purely so the
% figure can shade the nucleus and the two cytoplasm flanks in different colours.
left_search_start = max(1, nuc_col - 3*nuc_r_px);
left_search_end   = max(1, nuc_col - nuc_r_px - 1);
if left_search_start < left_search_end
    left_search_vals      = row_mean_img(left_search_start:left_search_end);
    [~, left_max_idx]     = max(left_search_vals);
    left_max_col          = left_search_start + left_max_idx - 1;
    mean_left_flank_start = max(1, left_max_col - nuc_r_px);
    mean_left_flank_end   = min(img_cols, left_max_col + nuc_r_px);
    mean_left_flank_end   = min(mean_left_flank_end, nuc_col - nuc_r_px - 1);
else
    mean_left_flank_start = left_search_start;
    mean_left_flank_end   = left_search_end;
end

right_search_start = min(img_cols, nuc_col + nuc_r_px + 1);
right_search_end   = min(img_cols, nuc_col + 3*nuc_r_px);
if right_search_start < right_search_end
    right_search_vals      = row_mean_img(right_search_start:right_search_end);
    [~, right_max_idx]     = max(right_search_vals);
    right_max_col          = right_search_start + right_max_idx - 1;
    mean_right_flank_start = max(1, right_max_col - nuc_r_px);
    mean_right_flank_end   = min(img_cols, right_max_col + nuc_r_px);
    mean_right_flank_start = max(mean_right_flank_start, nuc_col + nuc_r_px + 1);
else
    mean_right_flank_start = right_search_start;
    mean_right_flank_end   = right_search_end;
end

% --- y-limits: max from nucleus+flank, min from global ---
mean_plus_sem  = row_mean_img + row_sem_img;
mean_minus_sem = row_mean_img - row_sem_img;

% y-axis top is driven by the interesting region (flanks + nucleus) rather than
% the whole trace, so tall noise in the far periphery doesn't squash the plot;
% y-axis bottom uses the global minimum. Fall back to the full valid range if the
% flank region came out empty.
region_lo = max(1, round(mean_left_flank_start));
region_hi = min(img_cols, round(mean_right_flank_end));

region_plus_sem = mean_plus_sem(region_lo:region_hi);
region_plus_sem = region_plus_sem(~isnan(region_plus_sem));
if isempty(region_plus_sem); region_plus_sem = mean_plus_sem(valid_mean_img); end

sig_max = max(region_plus_sem);
sig_min = min(mean_minus_sem(valid_mean_img));
sig_pad = (sig_max - sig_min) * 0.1;

% --- Partition the linescan x-axis into 7 contiguous segments ---
% Ordered left->right: outer cyto | left flank | gap | NUCLEUS | gap | right
% flank | outer cyto. Each segment gets its own colour so the trace visually
% encodes which region each stretch belongs to. The boundaries are the flank
% windows and nucleus edges computed above.
x_ls = x_crop_img;

outer_left_idx  = x_ls <  mean_left_flank_start;
left_flank_idx  = x_ls >= mean_left_flank_start & x_ls <= mean_left_flank_end;
mid_left_idx    = x_ls >  mean_left_flank_end   & x_ls <  nuc_left;
nuc_idx         = x_ls >= nuc_left              & x_ls <= nuc_right;
mid_right_idx   = x_ls >  nuc_right             & x_ls <  mean_right_flank_start;
right_flank_idx = x_ls >= mean_right_flank_start & x_ls <= mean_right_flank_end;
outer_right_idx = x_ls >  mean_right_flank_end;

cyto_col  = [0.15 0.15 0.15];   % was white — invisible against a white figure background
flank_col = [0.0 0.75 0.6];     % slightly darkened teal for better contrast on white
nuc_col_c = [0.0 0.45 0.85];    % slightly darkened blue for better contrast on white

seg_indices = {outer_left_idx, left_flank_idx, mid_left_idx, nuc_idx, mid_right_idx, right_flank_idx, outer_right_idx};
seg_colors  = {cyto_col, flank_col, cyto_col, nuc_col_c, cyto_col, flank_col, cyto_col};

% =========================================================================
% Figure 1: Average projection + linescan trace panel
% =========================================================================
ax_sig_left   = 0.12;
ax_sig_width  = 0.78;

fig_height_img   = 150;
ax_sig_h_px      = 125;
slider_h_px      = 60;
top_margin_px    = 60;
fig_height_total = fig_height_img + ax_sig_h_px + slider_h_px + top_margin_px;

ax_img_aspect = img_rows / img_cols;
ax_img_w_px   = fig_height_img / ax_img_aspect;
fig_width     = ax_img_w_px / ax_sig_width;

% All positions normalized to total figure height
ax_img_bottom = slider_h_px    / fig_height_total;
ax_img_height = fig_height_img / fig_height_total;
ax_sig_height = ax_sig_h_px    / fig_height_total;
ax_sig_bottom = ax_img_bottom  + ax_img_height;

% Title sits inside ax_sig so top of ax_sig needs a little room
% Total used = ax_sig_bottom + ax_sig_height should be < 1
% = (slider + img + sig) / total = (60+300+250)/610 = 0.98 — fine

fig_main = figure('Color', 'w', 'Position', [100 100 round(fig_width) fig_height_total]);

ax_sig = axes('Parent', fig_main, ...
              'Position', [ax_sig_left ax_sig_bottom ax_sig_width ax_sig_height]);
ax_img = axes('Parent', fig_main, ...
              'Position', [ax_sig_left ax_img_bottom ax_sig_width ax_img_height]);

hold(ax_sig, 'on');

% SEM shading around mean
x_valid = x_crop_img(valid_mean_img);
fill(ax_sig, [x_valid fliplr(x_valid)], ...
     [row_mean_img(valid_mean_img) + row_sem_img(valid_mean_img), ...
      fliplr(row_mean_img(valid_mean_img) - row_sem_img(valid_mean_img))], ...
     [0.45 0.45 0.45], 'FaceAlpha', 0.30, 'EdgeColor', 'none');

% Colored mean line by segment
for s = 1:numel(seg_indices)
    idx = seg_indices{s};
    if ~any(idx); continue; end
    iidx  = find(idx);
    i1    = max(1, iidx(1)-1);
    i2    = min(numel(x_ls), iidx(end)+1);
    xs    = x_ls(i1:i2);
    ys    = row_mean_img(i1:i2);
    valid = ~isnan(ys);
    if ~any(valid); continue; end
    plot(ax_sig, xs(valid), ys(valid), '-', 'Color', seg_colors{s}, 'LineWidth', 2.5);
end

hold(ax_sig, 'off');

sig_ylim_lo = floor((sig_min) * 10) / 10;
n_steps     = max(1, ceil((sig_max - sig_ylim_lo)/0.2 - 1e-9));
sig_ylim_hi = sig_ylim_lo + n_steps*0.2;
sig_yticks  = sig_ylim_lo + (0:n_steps)*0.2;

set(ax_sig, 'Color', 'w', 'XColor', 'none', 'YColor', 'k', 'Box', 'off', ...
            'XLim', [0.5 img_cols+0.5], 'YLim', [sig_ylim_lo, sig_ylim_hi], ...
            'YTick', sig_yticks, ...
            'TickDir', 'out', 'TickLength', [0.01 0.01]);
ylabel(ax_sig, 'Intensity', 'Color', 'k');
title(ax_sig, 'Horizontal linescan — mean ± SEM', 'Color', 'k');

% Image

% skip dark opening colors until luminance clears a floor
% (min_lum), THEN ramp white into what's left -- so the display goes white ->
% Finding the skip point by luminance (rather than a fixed fraction of the
% colormap) adapts automatically to how slowly/quickly a given colormap darkens.
% min_lum   = 0.15;                                  % never show a color darker than this
% lum       = mean(cmap_thermal, 2);
% k         = find(lum >= min_lum, 1, 'first');
% if isempty(k); k = size(cmap_thermal, 1); end       % fallback: colormap never clears the floor
% cmap_trim = cmap_thermal(k:end, :);
% 
% white_frac = 0.18;
% n_ramp     = max(4, round(white_frac * size(cmap_trim, 1)));
% ramp       = [linspace(1, cmap_trim(1,1), n_ramp)', ...
%               linspace(1, cmap_trim(1,2), n_ramp)', ...
%               linspace(1, cmap_trim(1,3), n_ramp)'];
% cmap_black   = [ramp; cmap_trim];
clim_display = [0, clim(2)];

% Make dimmer pixels more transparent: alpha scales with each pixel's own
% intensity (normalized to the display range), so the low-intensity background
% near the mask edge fades smoothly into whatever's behind the axes instead of
% showing a hard/harsh boundary, while bright cell-body pixels stay opaque.
% alpha_gamma < 1 reaches full opacity sooner (only the genuinely dim pixels stay
% transparent); = 1 is a plain linear ramp; > 1 keeps more of the image faded.
alpha_gamma = 0.5;
alpha_mask  = min(1, max(0, avg_smooth / clim(2))) .^ alpha_gamma;

himg = imagesc('Parent', ax_img, 'CData', avg_smooth, clim_display);
set(himg, 'AlphaData', alpha_mask);
colormap(ax_img, 'hot');
set(ax_img, 'Color', 'w');   % shows through wherever alpha_mask < 1
axis(ax_img, 'off');
set(ax_img, 'DataAspectRatio', [1 1 1]);
set(ax_img, 'XLim', [0.5 img_cols+0.5], 'YLim', [0.5 img_rows+0.5]);
% title(ax_img, 'Average projection', 'Color', 'k');

hold(ax_img, 'on');
th = linspace(0, 2*pi, 100);
rc = rr_img;   % nucleus row in cropped coords (matches cc_img treatment for columns)
% plot(ax_img, cc_img + mean_nuc_r_px*cos(th), rc + mean_nuc_r_px*sin(th), ...
%      'w--', 'LineWidth', 2);
anno_color    = [0.15 0.15 0.15];
% xline(ax_img, cc_img, '--', 'Color',anno_color, 'LineWidth', 0.5);
yline(ax_img, nuc_row_img, ':', 'Color', anno_color, 'LineWidth', 0.5);
hold(ax_img, 'off');

linkaxes([ax_img ax_sig], 'x');

% Nucleus vertical annotation line spanning both panels, down to the bottom of
% the image panel (not just the top of it, where the trace panel ends).

x_frac        = (cc_img - 0.5) / img_cols;
line_x        = ax_sig_left + x_frac * ax_sig_width;
line_y_top    = ax_sig_bottom + ax_sig_height;
line_y_bottom = ax_img_bottom;
annotation(fig_main, 'line', [line_x line_x], [line_y_bottom line_y_top], ...
           'Color', anno_color, 'LineStyle', ':', 'LineWidth', 0.5);

% Sliders
slider_bottom = 0.005;
slider_h_norm = 0.03;

uicontrol('Parent', fig_main, 'Style', 'text', ...
          'Units', 'normalized', 'Position', [0.01 slider_bottom 0.08 slider_h_norm], ...
          'String', 'Brightness', 'BackgroundColor', 'w', 'ForegroundColor', 'k', 'FontSize', 8);
sl_bright = uicontrol('Parent', fig_main, 'Style', 'slider', ...
                       'Units', 'normalized', 'Position', [0.10 slider_bottom 0.18 slider_h_norm], ...
                       'Min', -1, 'Max', 1, 'Value', 0, 'Callback', @on_bright_slider, ...
                       'Tag', 'climBrightSlider');
sp_bright = uicontrol('Parent', fig_main, 'Style', 'edit', ...
                       'Units', 'normalized', 'Position', [0.29 slider_bottom 0.06 slider_h_norm], ...
                       'String', '0', 'BackgroundColor', [0.92 0.92 0.92], ...
                       'ForegroundColor', 'k', 'FontSize', 8, 'Callback', @on_bright_edit, ...
                       'Tag', 'climBrightEdit');

uicontrol('Parent', fig_main, 'Style', 'text', ...
          'Units', 'normalized', 'Position', [0.38 slider_bottom 0.08 slider_h_norm], ...
          'String', 'Contrast', 'BackgroundColor', 'w', 'ForegroundColor', 'k', 'FontSize', 8);
sl_contrast = uicontrol('Parent', fig_main, 'Style', 'slider', ...
                          'Units', 'normalized', 'Position', [0.47 slider_bottom 0.18 slider_h_norm], ...
                          'Min', 0.1, 'Max', 3, 'Value', 1, 'Callback', @on_contrast_slider, ...
                          'Tag', 'climContrastSlider');
sp_contrast = uicontrol('Parent', fig_main, 'Style', 'edit', ...
                          'Units', 'normalized', 'Position', [0.66 slider_bottom 0.06 slider_h_norm], ...
                          'String', '1', 'BackgroundColor', [0.92 0.92 0.92], ...
                          'ForegroundColor', 'k', 'FontSize', 8, 'Callback', @on_contrast_edit, ...
                          'Tag', 'climContrastEdit');

clim_range  = clim_display(2) - clim_display(1);
clim_center = mean(clim_display);

    function update_clim()
        brightness = get(sl_bright,   'Value');
        contrast   = get(sl_contrast, 'Value');
        new_range  = clim_range / contrast;
        new_center = clim_center + brightness * clim_range;
        set(ax_img, 'CLim', [new_center - new_range/2, new_center + new_range/2]);
    end

    function on_bright_slider(src, ~)
        set(sp_bright, 'String', sprintf('%.2f', src.Value));
        update_clim();
    end

    function on_bright_edit(src, ~)
        val = str2double(src.String);
        if ~isnan(val)
            val = max(-1, min(1, val));
            set(sl_bright, 'Value', val);
            set(src, 'String', sprintf('%.2f', val));
            update_clim();
        end
    end

    function on_contrast_slider(src, ~)
        set(sp_contrast, 'String', sprintf('%.2f', src.Value));
        update_clim();
    end

    function on_contrast_edit(src, ~)
        val = str2double(src.String);
        if ~isnan(val)
            val = max(0.1, min(3, val));
            set(sl_contrast, 'Value', val);
            set(src, 'String', sprintf('%.2f', val));
            update_clim();
        end
    end

% =========================================================================
% Figures 2-7 are optional (distributions / per-cell profiles). Skip them all
% when the caller only wants Figure 1.
% =========================================================================
if show_extra_figs

% =========================================================================
% Figure 2: All N/C ratios together with jitter
% =========================================================================
valid_r  = ~isnan(nc_ratios);
valid_m  = ~isnan(nc_ratios_mean);
valid_ls = ~isnan(col_nc_ratios);

jitter_amt = 0.12;
x_r   = 1 + (rand(sum(valid_r),  1) - 0.5) * jitter_amt;
x_m   = 2 + (rand(sum(valid_m),  1) - 0.5) * jitter_amt;
x_lsj = 3 + (rand(sum(valid_ls), 1) - 0.5) * jitter_amt;

figure; hold on;
plot(x_r,   nc_ratios(valid_r),      'o', 'MarkerFaceColor', [0.4 0.6 1],   'MarkerEdgeColor', 'none', 'MarkerSize', 8);
plot(x_m,   nc_ratios_mean(valid_m), 'o', 'MarkerFaceColor', [0.4 0.8 0.5], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
plot(x_lsj, col_nc_ratios(valid_ls), 'o', 'MarkerFaceColor', [0.9 0.5 0.2], 'MarkerEdgeColor', 'none', 'MarkerSize', 8);

plot([0.75 1.25], repmat(mean(nc_ratios(valid_r)),      1, 2), 'r-', 'LineWidth', 2);
plot([1.75 2.25], repmat(mean(nc_ratios_mean(valid_m)), 1, 2), 'r-', 'LineWidth', 2);
plot([2.75 3.25], repmat(mean(col_nc_ratios(valid_ls)), 1, 2), 'r-', 'LineWidth', 2);

yline(1, 'k:', 'LineWidth', 1.5);
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Radial', 'Simple mean', 'Linescan'}, 'XLim', [0.5 3.5]);
ylabel('N/C ratio');
title('N/C ratio comparison');
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 3: N/C ratio distribution (radial)
% =========================================================================
figure; hold on;
histogram(nc_ratios(valid_r), 'Normalization', 'probability', ...
          'FaceColor', [0.4 0.6 1], 'EdgeColor', 'none', 'NumBins', 10);
xline(1,                           'k:',  'LineWidth', 1.5);
xline(mean(nc_ratios(valid_r)),    'r-',  'LineWidth', 2);
xline(median(nc_ratios(valid_r)),  'w--', 'LineWidth', 1.5);
xlabel('N/C ratio (radial)'); ylabel('Probability');
title(sprintf('N/C ratio distribution — radial  (mean=%.2f, median=%.2f)', ...
              mean(nc_ratios(valid_r)), median(nc_ratios(valid_r))));
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 4: N/C ratio distribution (mean inside / mean outside)
% =========================================================================
figure; hold on;
histogram(nc_ratios_mean(valid_m), 'Normalization', 'probability', ...
          'FaceColor', [0.4 0.8 0.5], 'EdgeColor', 'none', 'NumBins', 10);
xline(1,                               'k:',  'LineWidth', 1.5);
xline(mean(nc_ratios_mean(valid_m)),   'r-',  'LineWidth', 2);
xline(median(nc_ratios_mean(valid_m)), 'w--', 'LineWidth', 1.5);
xlabel('N/C ratio (mean in / mean out)'); ylabel('Probability');
title(sprintf('N/C ratio distribution — simple mean  (mean=%.2f, median=%.2f)', ...
              mean(nc_ratios_mean(valid_m)), median(nc_ratios_mean(valid_m))));
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 5: Linescan N/C ratio distribution
% =========================================================================
figure; hold on;
histogram(col_nc_ratios(valid_ls), 'Normalization', 'probability', ...
          'FaceColor', [0.9 0.5 0.2], 'EdgeColor', 'none', 'NumBins', 10);
xline(1,                               'k:',  'LineWidth', 1.5);
xline(mean(col_nc_ratios(valid_ls)),   'r-',  'LineWidth', 2);
xline(median(col_nc_ratios(valid_ls)), 'w--', 'LineWidth', 1.5);
xlabel('N/C ratio (linescan)'); ylabel('Probability');
title(sprintf('Linescan N/C ratio distribution  (mean=%.2f, median=%.2f)', ...
              mean(col_nc_ratios(valid_ls)), median(col_nc_ratios(valid_ls))));
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 6: Individual column profiles — full canvas
% =========================================================================
figure('Color', 'k'); hold on;
for i = 1:N
    profile_i = col_profiles(:, i)';
    valid_i   = ~isnan(profile_i);
    if any(valid_i)
        plot(x_full(valid_i), profile_i(valid_i), ...
             'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
    end
end
plot(x_full(valid_col), col_mean_full(valid_col), 'w-', 'LineWidth', 2.5);
xline(cc_mid, 'w--', 'LineWidth', 1, 'Label', 'Nucleus');
xlabel('Canvas column index'); ylabel('Mean normalized intensity');
title('Column intensity profiles — individual cells + mean');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 7: Individual horizontal linescans — full canvas
% =========================================================================
figure('Color', 'k'); hold on;
for i = 1:N
    profile_i = row_profiles(:, i)';
    valid_i   = ~isnan(profile_i);
    if any(valid_i)
        plot(x_full(valid_i), profile_i(valid_i), ...
             'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
    end
end
fill([x_full(valid_row) fliplr(x_full(valid_row))], ...
     [row_mean_full(valid_row) + row_sem_full(valid_row), ...
      fliplr(row_mean_full(valid_row) - row_sem_full(valid_row))], ...
     [0.4 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(x_full(valid_row), row_mean_full(valid_row), 'w-', 'LineWidth', 2.5);
xline(cc_mid, 'w--', 'LineWidth', 1, 'Label', 'Nucleus');
xlabel('Canvas column index'); ylabel('Normalized intensity');
title('Horizontal linescan — individual cells + mean');
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
grid on;
cap_yticks(gca);

% =========================================================================
% Figure 8: Concentration metrics (nuclear fraction, Gini, entropy, and their
% linescan-restricted variants) as jittered points with group means. Those five
% live on a [0,1] scale and share the LEFT axis; r50 is a radius in pixels, so it
% gets its own RIGHT axis (last x position). Skipped if the struct predates them.
% =========================================================================
conc_data   = {nuc_sig_frac, gini_coef, entropy_coef, ls_gini_coef, ls_entropy_coef};
conc_labels = {'nuc frac', 'Gini', 'entropy', 'Gini (ls)', 'entropy (ls)'};
conc_have   = ~cellfun(@isempty, conc_data);
if any(conc_have) || ~isempty(radial_half_px)
    figure('Color', 'w'); hold on;
    jitter_amt = 0.14;
    n_conc     = numel(conc_data);
    r50_x      = n_conc + 1;                     % r50 sits just past the [0,1] metrics

    % --- Left axis: the [0,1] concentration metrics ---
    yyaxis left;
    for m = 1:n_conc
        if isempty(conc_data{m}); continue; end
        v  = conc_data{m}(:);
        v  = v(~isnan(v));
        if isempty(v); continue; end
        xm = m + (rand(numel(v),1) - 0.5) * jitter_amt;
        plot(xm, v, 'o', 'MarkerFaceColor', [0.30 0.45 0.75], ...
             'MarkerEdgeColor', 'none', 'MarkerSize', 7);
        plot([m-0.28 m+0.28], repmat(mean(v),1,2), '-', 'Color', [0.85 0.2 0.2], 'LineWidth', 2.5);
    end
    ylabel('metric value  (0 = spread, 1 = concentrated)');
    set(gca, 'YLim', [0 1], 'YColor', 'k');

    % --- Right axis: r50 in pixels ---
    all_labels = conc_labels;
    if ~isempty(radial_half_px)
        yyaxis right;
        hold on;
        v50 = radial_half_px(~isnan(radial_half_px));
        if ~isempty(v50)
            x50 = r50_x + (rand(numel(v50),1) - 0.5) * jitter_amt;
            plot(x50, v50, 'o', 'MarkerFaceColor', [0.90 0.55 0.10], ...
                 'MarkerEdgeColor', 'none', 'MarkerSize', 7);
            plot([r50_x-0.28 r50_x+0.28], repmat(mean(v50),1,2), '-', ...
                 'Color', [0.85 0.2 0.2], 'LineWidth', 2.5);
            ylabel('r50  (px)');
            set(gca, 'YColor', [0.75 0.45 0.05], 'YLim', [0 max(v50)*1.15]);
        end
        all_labels = [conc_labels, {'r50 (px)'}];
    end

    set(gca, 'XTick', 1:numel(all_labels), 'XTickLabel', all_labels, ...
             'XLim', [0.5 numel(all_labels)+0.5]);
    title('Concentration metrics — per cell + group mean');
    grid on;
end

% =========================================================================
% Figure 9: Radial intensity profile about the nucleus centre -- individual
% cells + population mean +/- SEM. Reads a peak-at-centre (nuclear-included) vs a
% dip-then-rise (nuclear-EXCLUDED) directly off the profile shape. Skipped if the
% struct predates the radial profile.
% =========================================================================
if ~isempty(radial_profiles)
    if ~isempty(radial_r_px)
        xr = radial_r_px(:)';
    else
        xr = (0:size(radial_profiles,1)-1) + 0.5;   % fallback if radius axis absent
    end
    rad_mean = mean(radial_profiles, 2, 'omitnan')';
    rad_sem  = (std(radial_profiles, 0, 2, 'omitnan') / sqrt(N))';
    valid_rr = ~isnan(rad_mean);

    figure('Color', 'k'); hold on;
    for i = 1:N
        pr = radial_profiles(:, i)';
        vi = ~isnan(pr);
        if any(vi)
            plot(xr(vi), pr(vi), 'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
            % Mark this cell's half-max radius ON its own curve: the point where
            % the profile first reaches half its own maximum. Sits at (r, max/2)
            % by construction, so the markers trace each cell's own half-max level.
            if ~isempty(radial_halfmax_px) && i <= numel(radial_halfmax_px) ...
                    && ~isnan(radial_halfmax_px(i))
                pv  = pr(vi);
                lvl = min(pv) + (max(pv) - min(pv))/2;   % matches RADIAL_HALF_MAX_RADIUS 'range'
                plot(radial_halfmax_px(i), lvl, 'o', ...
                     'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'w', ...
                     'MarkerSize', 6, 'LineWidth', 0.5);
            end
        end
    end
    if any(valid_rr)
        fill([xr(valid_rr) fliplr(xr(valid_rr))], ...
             [rad_mean(valid_rr) + rad_sem(valid_rr), ...
              fliplr(rad_mean(valid_rr) - rad_sem(valid_rr))], ...
             [0.4 0.6 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        plot(xr(valid_rr), rad_mean(valid_rr), 'w-', 'LineWidth', 2.5);
    end
    xlabel('Radius from nucleus centre (px)'); ylabel('Mean normalized intensity');
    if ~isempty(radial_halfmax_px) && any(~isnan(radial_halfmax_px))
        title(sprintf('Radial intensity profile — individual cells + mean  (o = half-max, median %.1f px)', ...
                      median(radial_halfmax_px(~isnan(radial_halfmax_px)))));
    else
        title('Radial intensity profile — individual cells + mean');
    end
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    grid on;
    cap_yticks(gca);
end

% =========================================================================
% Figure 10: Radial SUM profile -- total intensity per ring vs radius, individual
% cells + population mean. The sum rises with ring circumference, so its shape
% differs from the mean profile; it is the raw material for the cumulative below.
% =========================================================================
if ~isempty(radial_sums)
    if ~isempty(radial_r_px)
        xs = radial_r_px(:)';
    else
        xs = (0:size(radial_sums,1)-1) + 0.5;
    end
    rsum_mean = mean(radial_sums, 2, 'omitnan')';
    valid_rs  = ~isnan(rsum_mean);

    figure('Color', 'k'); hold on;
    for i = 1:N
        pr = radial_sums(:, i)';
        vi = ~isnan(pr);
        if any(vi)
            plot(xs(vi), pr(vi), 'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
        end
    end
    if any(valid_rs)
        plot(xs(valid_rs), rsum_mean(valid_rs), 'w-', 'LineWidth', 2.5);
    end
    xlabel('Radius from nucleus centre (px)'); ylabel('Total intensity in ring');
    title('Radial signal SUM per ring — individual cells + mean');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
    grid on;
    cap_yticks(gca);

    % =====================================================================
    % Figure 11: Cumulative radial distribution ("cumulative histogram of radial
    % sums"). For each cell, cumfrac(r) = cumsum(ring sums up to r) / total signal,
    % so every curve runs 0 -> 1 monotonically. The SHAPE is the discriminator:
    % nuclear-INCLUDED signal accumulates fast near the centre (curve bows up
    % early); nuclear-EXCLUDED signal accumulates late (curve stays low, then
    % climbs). Population mean overlaid.
    % =====================================================================
    cumfrac = NaN(size(radial_sums));
    for i = 1:N
        col = radial_sums(:, i);
        tot = sum(col, 'omitnan');
        if tot > 0
            col(isnan(col)) = 0;
            cumfrac(:, i) = cumsum(col) / tot;
        end
    end
    cf_mean  = mean(cumfrac, 2, 'omitnan')';
    valid_cf = ~isnan(cf_mean);

    figure('Color', 'k'); hold on;
    for i = 1:N
        cf = cumfrac(:, i)';
        vi = ~isnan(cf);
        if any(vi)
            plot(xs(vi), cf(vi), 'Color', [cmap(i,:) 0.4], 'LineWidth', 1);
        end
    end
    if any(valid_cf)
        plot(xs(valid_cf), cf_mean(valid_cf), 'w-', 'LineWidth', 2.5);
        % Mark r50: the radius where the MEAN cumulative curve crosses 0.5 (the
        % median radius of the signal), the scalar summary of this plot.
        kc = find(cf_mean >= 0.5, 1, 'first');
        if ~isempty(kc) && kc > 1
            c0 = cf_mean(kc-1); c1 = cf_mean(kc);
            if c1 > c0
                r50 = xs(kc-1) + (0.5 - c0)/(c1 - c0) * (xs(kc) - xs(kc-1));
            else
                r50 = xs(kc);
            end
            xline(r50, 'Color', [1 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5, ...
                  'Label', sprintf('r50 = %.1f px', r50), 'LabelVerticalAlignment', 'bottom');
        end
    end
    yline(0.5, 'w:', 'LineWidth', 1, 'Alpha', 0.5);
    xlabel('Radius from nucleus centre (px)');
    ylabel('Cumulative fraction of signal within radius');
    title('Cumulative radial signal — individual cells + mean');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'YLim', [0 1]);
    grid on;
    cap_yticks(gca);
end

end   % if show_extra_figs

end
% =========================================================================
function v = getfield_or_empty(s, name)
% Return s.(name) if present, else [] (so optional metrics/plots degrade
% gracefully for older structs that predate them).
    if isfield(s, name); v = s.(name); else; v = []; end
end
% =========================================================================
function cap_yticks(ax)
% Ensure the axis is visibly capped: force a tick at both the top and bottom
% y-limit, keeping whatever interior ticks are already there.
    yl = get(ax, 'YLim');
    yt = get(ax, 'YTick');
    yt = unique(round([yl(1), yt(yt > yl(1) & yt < yl(2)), yl(2)] * 1e6) / 1e6);
    set(ax, 'YTick', yt);
end
