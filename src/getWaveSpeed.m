function [vStat, tBins, xBins, vBins, Ebins] = getWaveSpeed(caseName, binWidth, tAvgWindow, xAvgWindow)
%GETWAVESPEED Bin a DARTESS run's wave-front trajectory and compute its speed.
%
%   [vStat, tBins, xBins, vBins, Ebins] = GETWAVESPEED(caseName, binWidth)
%   bins the front trajectory saved in <getTopResDir(caseName)>/drs.mat
%   (drs.xIwavpos, drs.t, drs.Ewavpos) and returns:
%     vStat - 1x3 summary [avgSpeedByTime, avgSpeedBySpace, peakSpeed], all
%             [m/s], over the (possibly windowed) analyzed trajectory:
%               avgSpeedByTime  = total distance / total time -- the usual
%                                 "average speed" (time-mean speed); equal to
%                                 the distance-weighted mean of the per-bin
%                                 speeds when weighted by TIME.
%                                 DEFAULT WINDOW (neither tAvgWindow nor
%                                 xAvgWindow supplied): anchored at the front's
%                                 own propagation ONSET, not at the run's t=0 --
%                                 see "Default analysis window" below. Pass an
%                                 explicit window to average over something
%                                 else (e.g. the full run including the
%                                 pre-propagation dead time: tAvgWindow =
%                                 [0, drs.t(end)]).
%               avgSpeedBySpace = per-bin speeds weighted by DISTANCE instead
%                                 of time, sum(v_i*dx_i)/sum(dx_i) (the
%                                 "space-mean speed" of kinematics/traffic-
%                                 flow theory) -- differs from avgSpeedByTime
%                                 whenever the speed varies along the
%                                 trajectory (it does here: fast initial
%                                 rise, then decline), since slow segments
%                                 contribute more bins' worth of TIME but
%                                 less DISTANCE, and vice versa for fast ones.
%               peakSpeed       = max(vBins), the single fastest bin.
%             Independent of binWidth's mode (time- or position-binning,
%             below) -- both feed the same underlying per-bin (v_i, dx_i,
%             dt_i) triples, just cut into bins along a different axis.
%             Caveat: avgSpeedBySpace is NaN (0/0) if the analyzed trajectory
%             covers zero net distance (e.g. a tAvgWindow that only catches
%             the front before it starts moving) -- correctly undefined,
%             since a distance-weighted average has no meaning over zero
%             distance; avgSpeedByTime stays well-defined (0) in that case.
%     tBins, xBins - the time [s] and position [m] of each bin's trailing
%             sample.
%     vBins - average front speed over each bin [m/s], =
%             -diff(position)/diff(time) between consecutive bin samples --
%             NOT an average of instantaneous per-step speeds; same
%             finite-difference-from-endpoints convention as plotWaveFront.m.
%     Ebins - MEAN front electric field [V/m] within each bin (all raw
%             per-step Ewavpos samples strictly between the previous bin edge
%             and this one, inclusive of this one) -- unlike xBins/tBins/
%             vBins, this is a true bin average, not a value sampled at the
%             trailing edge: the raw per-step |E| at the tracked front is
%             genuinely noisy (not just finely sampled), so decimating down
%             to the edge sample alone would carry that noise through rather
%             than smoothing it out. The first bin (no preceding edge) is
%             just its own single sample, same as xBins/tBins.
%
%   binWidth sets both the bin size AND, via comparison against the run's own
%   grid spacing drs.params.dx [m], which axis is binned:
%       binWidth <  dx : bin by TIME, width binWidth [s]      (e.g. 20e-6)
%       binWidth >= dx : bin by POSITION, width binWidth [m]  (e.g. 25)
%   This works because a physically sensible position-bin width can never be
%   finer than the grid (dx), while a time-bin width in seconds is always a
%   tiny number (~1e-9 to ~1e-3) next to dx (~O(1-10) m) -- so the comparison
%   unambiguously separates the two modes without a separate flag (contrast
%   plotWaveFront.m's explicit iBinByTime toggle).
%
%   [...] = GETWAVESPEED(caseName, binWidth, tAvgWindow) additionally
%   restricts the trajectory to samples with drs.t inside tAvgWindow, a 1x2
%   [tmin tmax] array [s] (order-independent), before binning. Pass [] to
%   skip this restriction.
%
%   [...] = GETWAVESPEED(caseName, binWidth, tAvgWindow, xAvgWindow) also
%   restricts to samples with drs.xIwavpos inside xAvgWindow, a 1x2
%   [xmin xmax] array [m] (order-independent). Both windows are optional and
%   combine as an intersection (AND) when both are given.
%
%   Default analysis window (neither tAvgWindow nor xAvgWindow supplied):
%   avgSpeedByTime is anchored at the front's own propagation ONSET -- the
%   first recorded sample where xIwavpos has actually left its starting
%   position (xIwavpos(1), normally x=0) -- through to the run's last sample,
%   NOT at the run's t=0. With a weak driving current the front can sit
%   stationary for a long time before it starts moving at all; averaging
%   total-distance/total-time from t=0 would silently dilute avgSpeedByTime
%   with that dead time. avgSpeedBySpace and peakSpeed need no such
%   adjustment -- dead-time bins already carry zero distance (dxSeg=0) and so
%   already contribute zero weight to avgSpeedBySpace, and a stationary bin
%   (speed=0) is never the max for peakSpeed. tBins/xBins/vBins/Ebins (the
%   full binned trajectory, e.g. for plotting) are NOT affected by this
%   default -- they still span the trajectory's full (possibly windowed)
%   extent, dead time included, exactly as before; only the avgSpeedByTime
%   scalar changes basis. Passing an EXPLICIT tAvgWindow/xAvgWindow (even an
%   empty one, []) always overrides this default and is used verbatim, same
%   as before -- e.g. tAvgWindow = [0, drs.t(end)] recovers the old
%   dead-time-inclusive average if that's what's wanted. If the front never
%   starts moving at all within the run, avgSpeedByTime falls back to the
%   full-trajectory calculation (correctly 0, not NaN or an error).
%
%   The first bin is split at the propagation onset (stationary sub-bin,
%   speed=0, followed by the propagating sub-bin), exactly as in
%   plotWaveFront.m, but ONLY when the (possibly windowed) trajectory being
%   binned still starts at the run's true t=0/x=0 -- the split is skipped
%   when tAvgWindow/xAvgWindow have cut off the true beginning, since there
%   is then no onset left inside the analyzed range to split at.
%
%   Shared by every script that reports a front speed, so that binning
%   and speed statistics mean the same thing everywhere.
%
%   Example:
%       [vStat, t, x, v, E] = getWaveSpeed('referenceCase', 20e-6);
%       plot(t*1e6, v*1e-6)   % speed [Mm/s] vs. time [us], time-binned
%       fprintf('avg (by time) = %.3g, avg (by space) = %.3g, peak = %.3g m/s\n', vStat)
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-08-11
% Last modified: 2026-10-09
%

narginchk(2,4)
if nargin < 3, tAvgWindow = []; end
if nargin < 4, xAvgWindow = []; end
noWindowGiven = isempty(tAvgWindow) && isempty(xAvgWindow);   % see "Default analysis window" above

%% Load
resdir  = getTopResDir(caseName);
matfile = [resdir 'drs.mat'];
if ~exist(matfile, 'file')
    error('getWaveSpeed:notFound', ...
        'No drs.mat found for caseName ''%s'' (expected %s).', caseName, matfile);
end
S   = load([resdir 'drs']);
drs = S.drs;

% Backward-compatible shim (files saved before 2026-08-09 used z-prefixed
% names for the horizontal coordinate; same convention as plotWaveFront.m).
if ~isfield(drs,'xIwavpos'), drs.xIwavpos = drs.zIwavpos; end
if ~isfield(drs,'xI'),       drs.xI       = drs.zI;       end

if ~isfield(drs.params,'dx')
    error('getWaveSpeed:noDx', ...
        'drs.params.dx missing -- cannot classify binWidth as time- or position-based.');
end
dx = drs.params.dx;

if ~isscalar(binWidth) || ~isreal(binWidth) || binWidth <= 0
    error('getWaveSpeed:badBinWidth', 'binWidth must be a positive real scalar.');
end
binByTime = binWidth < dx;

%% Full trajectory (before any windowing) -- kept to detect a true-start window below
mfrontFull = ~isnan(drs.xIwavpos);
xfrontFull = drs.xIwavpos(mfrontFull);
tfrontFull = drs.t(mfrontFull);
if isfield(drs,'Ewavpos')
    EfrontFull = drs.Ewavpos(mfrontFull);
else
    warning('getWaveSpeed:noEwavpos', 'drs.Ewavpos not found (older file?) -- Ebins will be NaN.');
    EfrontFull = nan(size(xfrontFull));
end

if numel(xfrontFull) < 2
    error('getWaveSpeed:tooFewSamples', ...
        'Fewer than 2 recorded front-position samples -- cannot compute a speed.');
end

% Propagation onset: first full-trajectory sample where the front has
% actually left its starting position (normally x=0) -- used below, when
% neither window was supplied, to anchor avgSpeedByTime past the
% pre-propagation dead time instead of at the run's t=0. Empty if the front
% never moves at all within the run (avgSpeedByTime then falls back to the
% full-trajectory 0/T = 0 result).
imoveFull = find(xfrontFull < xfrontFull(1), 1, 'first');

%% Apply optional windows (intersection)
keep = true(size(xfrontFull));
if ~isempty(tAvgWindow)
    if numel(tAvgWindow) ~= 2
        error('getWaveSpeed:badWindow', 'tAvgWindow must be a 1x2 array.');
    end
    keep = keep & (tfrontFull >= min(tAvgWindow)) & (tfrontFull <= max(tAvgWindow));
end
if ~isempty(xAvgWindow)
    if numel(xAvgWindow) ~= 2
        error('getWaveSpeed:badWindow', 'xAvgWindow must be a 1x2 array.');
    end
    keep = keep & (xfrontFull >= min(xAvgWindow)) & (xfrontFull <= max(xAvgWindow));
end

xfront = xfrontFull(keep);
tfront = tfrontFull(keep);
Efront = EfrontFull(keep);

if numel(xfront) < 2
    error('getWaveSpeed:emptyWindow', ...
        'Fewer than 2 front samples survive tAvgWindow/xAvgWindow -- window(s) too narrow or outside the run''s range.');
end

%% Bin edges (time- or position-based), same construction as plotWaveFront.m
if binByTime
    binKey   = tfront;
    binEdges = binKey(1):binWidth:binKey(end);
else
    binKey   = xfront;
    binEdges = binKey(1):-binWidth:binKey(end);
end
if numel(binEdges) < 2
    binEdges = [binKey(1) binKey(end)];   % binWidth wider than the whole (windowed) range: one bin
end
binEdges(end) = binKey(end);

Nedges = length(binEdges);
ibins  = zeros(1,Nedges);
for ib = 1:Nedges
    [~, ibins(ib)] = min(abs(binKey - binEdges(ib)));
end
ibins = unique(ibins, 'stable');   % guard against duplicate edges collapsing to the same sample (zero-width bin)

if numel(ibins) < 2
    error('getWaveSpeed:tooFewBins', ...
        'binWidth (with the given window) produces fewer than 2 distinct bin samples -- cannot compute a speed. Try a smaller binWidth or a wider window.');
end

xBins = xfront(ibins);
tBins = tfront(ibins);
Ebins = binMeanE(Efront, ibins);

%% Split the first bin at the propagation onset (stationary + propagating
% sub-bins), matching plotWaveFront.m -- only when this trajectory (after
% windowing) still starts at the run's true beginning, i.e. windowing hasn't
% already cut the onset out of the analyzed range.
startsAtTrueBeginning = (xfront(1) == xfrontFull(1)) && (tfront(1) == tfrontFull(1));
if startsAtTrueBeginning
    imove = find(xfront < xfront(1), 1, 'first');
    % Only insert a split point when the onset actually falls strictly
    % INSIDE the first natural bin (between ibins(1) and ibins(2)) -- if the
    % dead time is longer than binWidth, the onset can land at or past
    % ibins(2) instead, in which case the natural bins already separate dead
    % time from motion on their own and inserting imove-1 would produce a
    % non-monotonic (or duplicate) edge sequence: negative dtSeg/dxSeg
    % differences downstream (a real, previously silent bug -- the old
    % Efront(ibins) convention tolerated a duplicate/out-of-order index
    % without erroring, just producing a bogus bin; binMeanE's index-range
    % averaging does not).
    if ~isempty(imove) && imove - 1 > ibins(1) && imove - 1 < ibins(2)
        ibins = [ibins(1), imove-1, ibins(2:end)];
        xBins = xfront(ibins);
        tBins = tfront(ibins);
        Ebins = binMeanE(Efront, ibins);
    end
end

%% Average speed per bin: -diff(position)/diff(time) between consecutive bin
% samples (front moves to negative x) -- NOT an average of instantaneous
% per-step speeds, same convention as plotWaveFront.m. The first bin (no
% preceding edge to diff against) is assigned the same speed as the first
% real interval, for stairstep-style consumption (cf. plotWaveFront.m).
vBinsTry = -diff(xBins) ./ diff(tBins);
vBins    = [vBinsTry(1) vBinsTry];

%% Summary stats: two averages (time-weighted and space-weighted) plus peak.
% dxSeg/dtSeg are the true per-segment distance/duration (unpadded, unlike
% vBins/xBins/tBins themselves) -- the correct weights for both averages.
dxSeg = -diff(xBins);   % positive: front moves to negative x
dtSeg = diff(tBins);    % positive: time increases

if noWindowGiven && ~isempty(imoveFull)
    % No window was requested: anchor avgSpeedByTime at the propagation
    % onset (see "Default analysis window" above) rather than diluting it
    % with the pre-propagation dead time. Computed directly from the raw
    % (unbinned) onset/end samples -- equivalent to summing dxSeg/dtSeg over
    % just the post-onset bins, since the total telescopes regardless of
    % binning, but avoids having to re-locate the onset bin index here.
    % xBins/tBins/vBins/Ebins themselves are untouched: they still span the
    % run's full trajectory for plotting purposes.
    avgSpeedByTime = (xfrontFull(imoveFull) - xfrontFull(end)) / (tfrontFull(end) - tfrontFull(imoveFull));
else
    avgSpeedByTime = sum(dxSeg) / sum(dtSeg);   % = total distance/total time (telescopes regardless of binning)
end
avgSpeedBySpace = sum(vBinsTry .* dxSeg) / sum(dxSeg);   % distance-weighted; differs from avgSpeedByTime whenever speed isn't constant; already onset-invariant (dead-time bins carry dxSeg=0), no adjustment needed
peakSpeed       = max(vBins);

vStat = [avgSpeedByTime, avgSpeedBySpace, peakSpeed];

end


function Ebins = binMeanE(Efront, ibins)
%BINMEANE Mean of Efront over each bin defined by consecutive edge indices
%   ibins (assumed strictly increasing, as produced by the nearest-sample
%   search above on a monotonic binKey). Bin i (i>=2) averages every raw
%   sample strictly after the previous edge through this edge, inclusive:
%   Efront(ibins(i-1)+1 : ibins(i)). The first bin has no preceding edge, so
%   it is just its own single sample -- same convention as xBins/tBins,
%   which are always edge samples, never bin means.
Ebins = zeros(size(ibins));
Ebins(1) = Efront(ibins(1));
for i = 2:numel(ibins)
    Ebins(i) = mean(Efront(ibins(i-1)+1 : ibins(i)));
end
end
