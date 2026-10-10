% DARTESS : Plot dart leader wave front position, speed and electric field
% for one or more saved runs (simCases), overlaid. 2x2 layout: (1) front
% position vs. time, (2) front speed vs. position, (3) front speed vs. time,
% (4) |E| at the front vs. position. Panels 2 and 4 mark every loaded case's
% own interior-junction positions as light dashed reference lines.
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-06-06
% Last modified: 2026-10-09
%
clc; clear all

% getTopResDir.m/getWaveSpeed.m live in ./src, relative to the folder this
% script is run from.
if exist('src','dir'), addpath('src'); end

%% Cases to load and overlay
% caseName of each saved run to overlay (folder names under getTopResDir.m's
% root). The default is dartess.m's own default run.
simCases = {'referenceCase'};

Nsims = length(simCases);

%% Style
tunit = 1e6;  % Convert to us
vunit = 1e-6; % Convert to Mm/s
Eunit = 1e-5; % Convert to 1e5 V/m
LFS   = 24;
FS    = 20;
LW0    = 7;    % first curve's line width -- a bit thick to start
LWstep = 0.85; % each subsequent curve 15% thinner than the last, so an
               % overlapping earlier (thicker) curve still peeks out from
               % under a later (thinner) one instead of being fully hidden
lineWidths = LW0 * LWstep.^(0:Nsims-1);
cmap  = lines(Nsims); % distinct color per simulation

% Wave front speed estimate: bin the front trajectory either by position or
% by time before differencing (the "wide bins" method of TERS's plot_v_vs_z.m),
% then take diff(position)/diff(time) over each bin. The
% same bins are used to smooth the electric field at the wave front.
% Binning/speed calculation itself delegated to getWaveSpeed.m -- a single
% binWidth is all it needs: getWaveSpeed picks time- vs. position-binning by
% comparing binWidth against the run's own dx (binWidth < dx -> time [s];
% binWidth >= dx -> position [m]), so set binWidth in whichever unit is
% wanted directly, e.g. 20e-6 for 20us-wide time bins or 25 for 25m-wide
% position bins. Use narrower bins to resolve a short-lived transient.
binWidth = 20e-6;

f3 = figure(3); clf
f3.Position = [100 100 1100 850];

%% Pass 1: load every case once and build its legend label
% Loaded here (not re-loaded per case further down) so panel 1's front
% position vs. time -- the one panel that plots drs fields directly rather
% than through getWaveSpeed.m -- has what it needs without a second load.
% Labeling is its own complete pass, before any plotting or console output,
% so that the collision fallbacks below can rewrite every label before
% anything is printed or drawn.
drsAll     = cell(1,Nsims);
simCasestr = cell(1,Nsims);
ItipAll    = nan(1,Nsims);
xmarkAll   = cell(1,Nsims);
EambStr    = cell(1,Nsims);   % per-case ambient-field label, for the second-level fallback below
XminStr    = cell(1,Nsims);   % per-case domain label, for the fallback below -- the DOMAIN is a
                              % first-class property in this figure, since it sets where the
                              % EDGE early-stop cuts each front off (dibuf = 5% of the domain)

for j = 1:Nsims
    S = load([getTopResDir(simCases{j}) 'drs']);
    drs = S.drs;
    % Backward compatibility: files saved before 2026-08-09 used z-prefixed
    % names for the horizontal coordinate (renamed to x).
    if ~isfield(drs,'xIwavpos'), drs.xIwavpos = drs.zIwavpos; end
    if ~isfield(drs,'xI'), drs.xI = drs.zI; end
    if ~isfield(drs.params,'xjunc')
        if isfield(drs.params,'zjunc')
            drs.params.xjunc = drs.params.zjunc;
        elseif isfield(drs.params,'zbrnch')
            drs.params.xjunc = drs.params.zbrnch;
        end
    end
    drsAll{j} = drs;

    % Tip current: junction entry 1 since V42; older files carry params.Ilead
    if isfield(drs.params,'Ijunc')
        Itip = drs.params.Ijunc(1);
    else
        Itip = drs.params.Ilead;
    end
    ItipAll(j) = Itip;

    % Ambient field, Eamb(x) = E0 + s*x, read through src/getEambE0s.m so every saved
    % convention labels the same way (v0.102+ [E0 s]; v0.97-v0.101 [Val_at_xmax, Val_at_xmin];
    % older files have no ambient field).
    [eaE0, eaS] = getEambE0s(drs);
    if eaS == 0
        EambStr{j} = ['E_{amb} = ' sprintf('%.3g', eaE0*1e-3) ' kV/m'];
    else
        EambStr{j} = ['E_0 = ' sprintf('%.3g', eaE0*1e-3) ' kV/m, s = ' sprintf('%.3g', eaS) ' V/m^2'];
    end

    % Legend label: when all loaded cases share the same Ilead (as with a Scenario I/II/III
    % comparison), an Ilead-only label is useless -- every curve would read identically. Label
    % by scenario instead, inferred from the junction count (2/3/4 -> Scenario I/II/III, the
    % structural signature of the three example scenarios in dartess.m); falls
    % back to the Ilead label for anything that doesn't match.
    if isfield(drs.params,'xjunc')
        switch numel(drs.params.xjunc)
            case 2
                simCasestr{j} = 'Scenario I';
            case 3
                simCasestr{j} = 'Scenario II';
            case 4
                simCasestr{j} = 'Scenario III';
            otherwise
                simCasestr{j} = ['I_{lead} = ' sprintf('%g',Itip) ' A'];
        end
    else
        simCasestr{j} = ['I_{lead} = ' sprintf('%g',Itip) ' A'];
    end

    XminStr{j} = ['x_{min} = ' sprintf('%g', min(drs.xI)) ' m'];

    % Junction positions, for the reference lines on the vs-position panels.
    % Entry 1 is the leader tip at xmax (the domain edge), so it is skipped --
    % a line there would just be the axis. Kept per case (not overwritten) so
    % that comparisons across structurally different cases (e.g. Scenario
    % I/II/III, which have different junction layouts) mark every case's own
    % junctions, not just the last one loaded.
    if isfield(drs.params,'xjunc')
        xmarkAll{j} = drs.params.xjunc(2:end);
    else
        xmarkAll{j} = [];   % no junction info at all (pre-V41 file)
    end
end

% Inverse collision, e.g. an Ilead SWEEP where every leg happens to be
% Scenario I: the scenario-based label above comes out identical for every
% case, just as useless as the Ilead-only label was in the case it was
% designed to fix. Fall back to Ilead-based labels for all cases instead.
if Nsims > 1 && numel(unique(simCasestr)) == 1
    for j = 1:Nsims
        simCasestr{j} = ['I_{lead} = ' sprintf('%g',ItipAll(j)) ' A'];
    end
end

% Second-level collision: an ambient-field
% sweep is structurally ONE scenario at ONE Ilead, so the scenario label and the
% Ilead fallback above both come out identical for every leg -- two useless labels
% in a row. Fall back to Eamb, but only when it actually varies across the loaded
% set, so this never fires on a genuine Ilead sweep (where every leg has Eamb = 0
% and this would just replace a useful label with an identical one).
if Nsims > 1 && numel(unique(simCasestr)) == 1 && numel(unique(EambStr)) > 1
    simCasestr = EambStr;
end

% Third-level collision: a MIXED set in which several attributes vary at once
% -- e.g. a current-driven run loaded against ambient-field-driven runs at
% different currents AND different domains. Neither fallback above handles it:
% the Ilead fallback leaves same-current legs reading identically, and the Eamb
% fallback never fires because it requires the labels to be ALL identical.
%
% Rule: append every attribute that actually VARIES across the loaded set, in
% order. Two notes on why it is written this way rather than as another
% "fall back only if still ambiguous" step:
%   * It appends to EVERY label, not just the colliding ones, so the legend
%     never mixes two label formats.
%   * x_min is appended whenever it varies even if the labels are already
%     unique: two labels can be distinct yet nearly indistinguishable at a
%     glance while the runs behind them differ by kilometers of domain, and in
%     a figure whose panels plot front POSITION a reader must not have to infer
%     a different domain.
if Nsims > 1 && numel(unique(simCasestr)) < Nsims
    varyAttr = {EambStr, XminStr};
    for a = 1:numel(varyAttr)
        if numel(unique(varyAttr{a})) > 1
            for j = 1:Nsims
                simCasestr{j} = [simCasestr{j} ', ' varyAttr{a}{j}];
            end
        end
    end
end

%% Pass 2: bin each case's front trajectory, report, and plot
pX = gobjects(1,Nsims);
for j=1:Nsims
    drs = drsAll{j};

    %% Bin wave front trajectory to obtain smoother speed and field estimates
    % (delegated to getWaveSpeed.m)
    [vStat, tbins, xbins, vbins, Ebins] = getWaveSpeed(simCases{j}, binWidth);

    % Report max and average front speed
    maxSpeed = vStat(3)*vunit;   % peakSpeed
    avgSpeed = vStat(1)*vunit;   % avgSpeedByTime -- same total-distance/total-time convention as before
    disp([simCasestr{j} ': max speed = ' sprintf('%1.3f',maxSpeed) ' Mm/s, average speed = ' sprintf('%1.3f',avgSpeed) ' Mm/s']);

    LWj = lineWidths(j);


    %% Plot
    subplot(2,2,1) % Wave front position versus time
    pX(j) = plot(drs.t*tunit, drs.xIwavpos, '-', 'Color', cmap(j,:), 'LineWidth', LWj);
    if j==1
        hold on
    end
    if j==Nsims
        hold off
        set(gca,'YDir','reverse','FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        xlabel('Time (\mus)','FontSize',LFS)
        ylabel('Wave front position (m)','FontSize',LFS)
        % The legend inherits the axes FontSize (FS = 20), which is fine for a
        % few short labels but covers most of panel 1 once many long labels are
        % loaded. Scale it down by BOTH drivers of legend area -- the number of
        % entries and the longest label -- with a floor so it stays readable.
        % A 2-case set is unaffected (both factors clamp at 1).
        legFS = FS * min(1, 5/Nsims) * min(1, 38/max(cellfun(@numel, simCasestr)));
        legend(pX, simCasestr, 'Location', 'Best', 'FontSize', max(9, round(legFS)))
    end

    subplot(2,2,2) % Wave front speed versus position
    stairs(xbins, vbins*vunit, '-', 'Color', cmap(j,:), 'LineWidth', LWj);
    if j==1
        hold on
    end
    if j==Nsims
        for jj = 1:Nsims   % mark every loaded case's own junction positions
            for ib = 1:numel(xmarkAll{jj})
                line([xmarkAll{jj}(ib) xmarkAll{jj}(ib)], ylim, 'Color',[0.6 0.6 0.6], ...
                     'LineStyle','--', 'HandleVisibility','off')
            end
        end
        hold off
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on','XDir','reverse')
        xlabel('Wave front position (m)','FontSize',LFS)
        ylabel('Wave front speed (Mm/s)','FontSize',LFS)
    end

    subplot(2,2,3) % Wave front speed versus time
    stairs(tbins*tunit, vbins*vunit, '-', 'Color', cmap(j,:), 'LineWidth', LWj);
    if j==1
        hold on
    end
    if j==Nsims
        hold off
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
        xlabel('Time (\mus)','FontSize',LFS)
        ylabel('Wave front speed (Mm/s)','FontSize',LFS)
    end

    subplot(2,2,4) % Electric field at the wave front versus position
    % Ebins is the MEAN |E| within each bin (getWaveSpeed.m), not a raw
    % instantaneous sample -- the full-resolution drs.Ewavpos trace is
    % genuinely noisy step-to-step, so some form of averaging is needed for
    % a legible plot.
    stairs(xbins, Ebins*Eunit, '-', 'Color', cmap(j,:), 'LineWidth', LWj);
    if j==1
        hold on
    end
    if j==Nsims
        for jj = 1:Nsims   % mark every loaded case's own junction positions
            for ib = 1:numel(xmarkAll{jj})
                line([xmarkAll{jj}(ib) xmarkAll{jj}(ib)], ylim, 'Color',[0.6 0.6 0.6], ...
                     'LineStyle','--', 'HandleVisibility','off')
            end
        end
        hold off
        set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on','XDir','reverse')
        xlabel('Wave front position (m)','FontSize',LFS)
        ylabel('|E| at wave front (kV/cm)','FontSize',LFS)
    end

end
