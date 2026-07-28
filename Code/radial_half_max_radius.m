function [r, info] = radial_half_max_radius(prof, r_px, opts)
% RADIAL_HALF_MAX_RADIUS  Radius at which a radial profile first reaches its
% half-maximum level, scanning outward from the nucleus centre.
%
% Shared by AVERAGE_CELL_INTENSITY (which stores it as results.radial_halfmax_px),
% BACKFILL_RADIAL_METRICS (which adds it to already-computed results) and
% PLOT_RADIAL_BY_CONDITION (which marks it on the curves), so the three can never
% disagree about what "half max" means.
%
% The calculation is deliberately plain: take the profile as it is, find its
% maximum, and return the first radius at which it crosses half that value.
% No ring trimming, no smoothing -- every valid ring is used exactly as measured.
%
% Inputs:
%   prof : radial profile (a column of results.radial_profiles). NaN rings --
%          radii the cell does not reach -- are ignored.
%   r_px : radius at each ring centre (results.radial_r_px).
%   opts : optional struct (or, for back-compatibility, the level_mode string):
%            .level_mode  'max' (default) = half the maximum, peak/2
%                         'range'         = baseline + (peak-baseline)/2, with the
%                                           baseline taken as the lo_pct percentile
%            .smooth      moving-average window in rings, odd, 1 = OFF (default).
%                         WARNING: smoothing blurs a SHARP profile proportionally
%                         more than a broad one, so it compresses real differences
%                         between conditions that differ in sharpness. Measured on
%                         simulated profiles with half-max radii of ~2.4 vs 3.5 px,
%                         a 3-ring window recovered only 79% of the true difference
%                         against 91% unsmoothed. Left available, but off.
%            .lo_pct      baseline percentile for 'range' [10]; ignored by 'max'.
%
% Direction-agnostic: the scan starts at r = 0 and returns the FIRST crossing, so
% for a centre-peaked (nuclear-included) profile that is the crossing on the way
% DOWN, and for a centre-dipped (nuclear-excluded) one it is the crossing on the
% way UP. One rule, no special-casing.
%
% NOTE on 'max': if a profile sits on a pedestal higher than half its own peak, the
% level falls below the whole curve, nothing crosses it, and the result is NaN.
% That is a property of the level choice, not an error -- check
% sum(isnan(results.radial_halfmax_px)) if a batch comes back unexpectedly empty.
%
% Outputs:
%   r    : the interpolated radius (same units as r_px), or NaN if the profile has
%          fewer than two valid rings or never crosses the level.
%   info : struct with .level, .baseline, .peak, .n_rings_used and .profile_used --
%          for diagnosing a suspicious value.

    if nargin < 3 || isempty(opts); opts = struct(); end
    if ~isstruct(opts)
        opts = struct('level_mode', char(opts));      % back-compat: 3rd arg was a string
    end
    if ~isfield(opts, 'level_mode'); opts.level_mode = 'max'; end
    if ~isfield(opts, 'smooth');     opts.smooth     = 1;     end
    if ~isfield(opts, 'lo_pct');     opts.lo_pct     = 10;    end

    info = struct('level', NaN, 'baseline', NaN, 'peak', NaN, ...
                  'n_rings_used', 0, 'profile_used', []);

    prof = prof(:);
    r_px = r_px(:);

    % Use every ring the cell actually reaches. Nothing else is discarded.
    keep = ~isnan(prof);
    prof = prof(keep);
    r_px = r_px(keep);
    info.n_rings_used = numel(prof);
    if numel(prof) < 2
        r = NaN;
        return;
    end

    % Optional smoothing (off by default -- see the warning above).
    prof_s = movmean_nan(prof, opts.smooth);
    info.profile_used = prof_s;

    peak = max(prof_s);
    switch lower(opts.level_mode)
        case 'max'
            info.peak     = peak;
            info.baseline = 0;
            if ~(peak > 0); r = NaN; return; end
            level = peak / 2;
        case 'range'
            base = pctl_local(prof_s, opts.lo_pct);
            info.peak     = peak;
            info.baseline = base;
            if ~(peak > base); r = NaN; return; end     % flat profile: no half-max
            level = base + (peak - base) / 2;
        otherwise
            error('radial_half_max_radius: level_mode must be ''max'' or ''range''.');
    end
    info.level = level;

    % First crossing of the level, scanning outward, linearly interpolated.
    above = prof_s >= level;
    k     = find(above ~= above(1), 1, 'first');
    if isempty(k)
        r = NaN;
        return;
    end

    p0 = prof_s(k-1); p1 = prof_s(k);
    if p1 ~= p0
        r = r_px(k-1) + (level - p0) / (p1 - p0) * (r_px(k) - r_px(k-1));
    else
        r = r_px(k);
    end
end

% =========================================================================
function y = movmean_nan(x, w)
% Centred moving average over w samples, ignoring NaN. w <= 1 returns x unchanged.
    if w <= 1
        y = x;
        return;
    end
    n = numel(x);
    h = floor(w / 2);
    y = NaN(n, 1);
    for i = 1:n
        lo = max(1, i - h);
        hi = min(n, i + h);
        v  = x(lo:hi);
        v  = v(~isnan(v));
        if ~isempty(v); y(i) = mean(v); end
    end
end

% =========================================================================
function q = pctl_local(x, p)
% Percentile without the Statistics Toolbox (linear interpolation on sorted data).
    x = x(:); x = x(~isnan(x));
    if isempty(x); q = NaN; return; end
    xs = sort(x); n = numel(xs);
    if n == 1; q = xs; return; end
    pos = (p/100) * n + 0.5;
    pos = min(max(pos, 1), n);
    lo  = floor(pos); hi = ceil(pos);
    if lo == hi
        q = xs(lo);
    else
        q = xs(lo) + (pos - lo) * (xs(hi) - xs(lo));
    end
end
