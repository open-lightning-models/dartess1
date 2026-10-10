% gen1Dmov generates a movie (or on-screen animation) overlaying one or more
% DARTESS dart leader simulations, one color per simulation (legName). The
% panel layout is (Q, I) / (E, rQ) for a constant-R run, or (Q, I) / (E, rQ)
% / (T, S) when the thermodynamic closure was active (iupdateR>=1). Panel 4
% shows the corona sheath radius rQ(x,t), read from the saved drs.rQsave.
% This works best when tsave is defined in dartess.m like this,
% for example: Ntsave = 101; tsave = linspace(0,tmax,Ntsave);
% where Ntsave will be the number of frames in the movie
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-06-06
% Last modified: 2026-10-09
%
clc; clear all

% getTopResDir.m lives in ./src, relative to the folder this script is run
% from.
if exist('src','dir'), addpath('src'); end

%% Input
% caseName of each saved run to overlay (folder names under getTopResDir.m's
% root). The default is dartess.m's own default run.
simCases = {'referenceCase'};
legName  = {'Reference case'}; % Legend names, one per simCases entry

iSaveMov = 0; % =1 to save movie
% movName: output movie file name (no extension). Optional -- if undefined,
% defaults to the name of the first input simulation (simCases{1}).
% movName = 'myMovie';


%% Load every simulation ONCE, up front
% Every frame below uses only what is loaded here, so each file is read once.
Nsims = length(simCases);
drsAll   = cell(1,Nsims);
EsaveAll = cell(1,Nsims); % |E| is not saved by dartess.m -- recomputed once per sim, not once per frame
xEAll    = cell(1,Nsims);
for j = 1:Nsims
    S   = load([getTopResDir(simCases{j}) 'drs']);
    drs = S.drs;
    % Backward compatibility: files saved before 2026-08-09 used z-prefixed
    % names for the horizontal coordinate (renamed to x). Patching the loaded drs/params here lets
    % everything below use only the new x-names, regardless of file era.
    if ~isfield(drs,'xI'), drs.xI = drs.zI; drs.xV = drs.zV; end
    if ~isfield(drs.params,'dx'), drs.params.dx = drs.params.dz; end
    if ~isfield(drs.params,'xjunc')
        if isfield(drs.params,'zjunc')
            drs.params.xjunc = drs.params.zjunc;
        elseif isfield(drs.params,'zbrnch')
            drs.params.xjunc = drs.params.zbrnch;
        end
    end
    drsAll{j}   = drs;
    xEAll{j}    = drs.xI(2:end-1);
end

% Every case is indexed by the SAME frame number ksave below (drs.Qsave(ksave,:),
% etc.), so every loaded case must share one tsave grid; otherwise a case would
% either index past the end of its own snapshots or silently pair up frames
% that are not at the same real time. Two runs meant for a side-by-side
% comparison often do NOT share one tsave grid, even with the same nominal
% Ntsave/tmax: an early stop at the domain edge trims tsave to the snapshots
% actually reached, and two runs can reach it at different times (e.g. a TE run
% and a QES run of the same configuration).
%
% So every OTHER simulation's snapshots are interpolated linearly in time onto
% simCases{1}'s tsave grid, clipped to whichever loaded run ended soonest
% (tCommonMax). This is exact wherever simCases{1}'s snapshot times are also
% snapshot times of the other run, and otherwise a linear interpolation between
% that run's own snapshots.
%
% Spatial grids are NOT interpolated -- an overlay only makes sense between
% simulations run on one common dx/domain in the first place (a side-by-side
% comparison at two different dx or domain sizes is not a fair overlay), so
% that case is rejected explicitly rather than silently misaligned.
dxAll   = cellfun(@(d) d.params.dx, drsAll);
xminAll = cellfun(@(d) min(d.xI),   drsAll);
xmaxAll = cellfun(@(d) max(d.xI),   drsAll);
if any(dxAll ~= dxAll(1)) || any(xminAll ~= xminAll(1)) || any(xmaxAll ~= xmaxAll(1))
    error('gen1Dmov:mismatchedGrid', ...
        ['Loaded cases do not share one spatial grid (dx = %s, xmin = %s, xmax = %s) -- ' ...
         'an overlay requires every simCases entry to share dx and domain.'], ...
        mat2str(dxAll), mat2str(xminAll), mat2str(xmaxAll));
end

tsaveAll = cellfun(@(d) d.tsave, drsAll, 'UniformOutput', false);
tsaveNumelAll = cellfun(@numel, tsaveAll);
if Nsims > 1 && (any(tsaveNumelAll ~= tsaveNumelAll(1)) || ~isequal(tsaveAll{:}))
    tCommonMax = min(cellfun(@(t) max(t), tsaveAll));
    tsaveRef   = drsAll{1}.tsave;
    tsaveRef   = tsaveRef(tsaveRef <= tCommonMax);
    fprintf(['gen1Dmov: loaded cases do not share one tsave grid -- interpolating onto\n' ...
             'simCases{1}''s own grid, clipped to %.4g us (the shortest loaded run).\n'], tCommonMax*1e6);
    for j = 1:Nsims
        drs = drsAll{j};
        if isequal(drs.tsave, tsaveRef), continue; end
        drs.Qsave  = interp1(drs.tsave, drs.Qsave,  tsaveRef);
        drs.Isave  = interp1(drs.tsave, drs.Isave,  tsaveRef);
        drs.Vsave  = interp1(drs.tsave, drs.Vsave,  tsaveRef);
        drs.rQsave = interp1(drs.tsave, drs.rQsave, tsaveRef);
        if isfield(drs,'Tsave')
            drs.Tsave  = interp1(drs.tsave, drs.Tsave,  tsaveRef);
            drs.Ssave  = interp1(drs.tsave, drs.Ssave,  tsaveRef);
            drs.Eksave = interp1(drs.tsave, drs.Eksave, tsaveRef);
        end
        drs.tsave = tsaveRef;
        drsAll{j} = drs;
    end
end

% |E| is not saved by dartess.m -- recomputed once per sim (from the final,
% possibly time-interpolated, Vsave), not once per frame.
for j = 1:Nsims
    EsaveAll{j} = -diff(drsAll{j}.Vsave, 1, 2)/drsAll{j}.params.dx;
end

%% Grid and domain, taken from the first simulation
% A frame-by-frame overlay only makes sense when every simCases entry shares
% one domain and one tsave grid, so the first simulation's grid/timebase is
% used to drive the whole movie (matching the original design).
drs1     = drsAll{1};
Ntsave   = length(drs1.tsave);
tsave    = drs1.tsave;
xI       = drs1.xI;
xmin     = min(xI);
xmax     = max(xI);
ldom     = xmax - xmin;
v        = drs1.params.v;
dx       = drs1.params.dx;
iupdateR = drs1.flags.iupdateR;

% Junction positions of the FIRST simulation, drawn as reference lines on every
% panel. Entry 1 is the leader tip at xmax (the domain edge), so it is skipped.
if isfield(drs1.params,'xjunc')
    xmark = drs1.params.xjunc(2:end);
else
    xmark = [];   % no junction info at all (pre-V41 file)
end
% line() rather than plot(): it adds to the axes regardless of the hold state,
% so it can be called after the panel's hold off / xlim / ylim block.
drawJunc = @() arrayfun(@(xb) line([xb xb], ylim, 'Color',[0.6 0.6 0.6], ...
                    'LineStyle','--', 'HandleVisibility','off'), xmark);

%% Axis limits, spanning ALL loaded simulations
% Taken across every loaded simulation (running min/max over j), so that no
% curve of an overlay is clipped by another simulation's range.
Qmin = Inf; Qmax = -Inf;
rQmin = Inf; rQmax = -Inf;
Imin = Inf; Imax = -Inf;
Eabs_all = [];
if iupdateR >= 1
    Tmin = Inf; Tmax = -Inf;
    Sall = [];
end
for j = 1:Nsims
    drs = drsAll{j};
    Qmin  = min(Qmin,  min(drs.Qsave(:)));  Qmax  = max(Qmax,  max(drs.Qsave(:)));
    rQmin = min(rQmin, min(drs.rQsave(:))); rQmax = max(rQmax, max(drs.rQsave(:)));
    Imin  = min(Imin,  drs.Ilims(1));       Imax  = max(Imax,  drs.Ilims(2));
    Eabs_all = [Eabs_all; abs(EsaveAll{j}(:))]; %#ok<AGROW> -- Nsims is small (a handful of overlaid runs)
    if iupdateR >= 1
        Tmin = min(Tmin, min(drs.Tsave(:))); Tmax = max(Tmax, max(drs.Tsave(:)));
        Sall = [Sall; drs.Ssave(:)]; %#ok<AGROW>
    end
end

Ilims = [Imin Imax];

dQp  = 0.05*(Qmax-Qmin);
Qlims = [Qmin-dQp Qmax+dQp];

% Corona sheath radius limits (panel 4), from the saved drs.rQsave.
if rQmax > rQmin
    drQp = 0.05*(rQmax-rQmin);
else
    drQp = max(0.1*rQmax, 1e-3);   % static-radius runs (rQ uniform everywhere): small fixed margin
end
rQlims = [rQmin-drQp rQmax+drQp];

% Electric field limits (|E| is not saved by dartess.m; it is computed above from Vsave)
Emax_val = max(Eabs_all);
dE_orders_of_magnitude = 5;
if Emax_val == 0; Emax_val = 1; end
Esig = Eabs_all(Eabs_all > (10^(dE_orders_of_magnitude-1))*Emax_val);
if isempty(Esig)
    Elogmin = ceil(log10(Emax_val)) - dE_orders_of_magnitude;
else
    Elogmin = floor(log10(min(Esig)));
end
Elogmax = ceil(log10(Emax_val));
Elims   = 10.^[Elogmin, Elogmax];

if iupdateR >= 1
    dTp  = max(0.05*(Tmax-Tmin), 100);
    Tlims = [max(Tmin-dTp,0) Tmax+dTp];

    dS_orders_of_magnitude = 8;
    Smax_val = max(Sall); if Smax_val == 0; Smax_val = 1; end
    Ssig = Sall(Sall > (10^(dS_orders_of_magnitude-1))*Smax_val);
    if isempty(Ssig)
        Slogmin = ceil(log10(Smax_val)) - dS_orders_of_magnitude;
    else
        Slogmin = floor(log10(min(Ssig(Ssig>0))));
    end
    Slogmax = ceil(log10(Smax_val));
    Slims   = 10.^[Slogmin, Slogmax];
end


%% Movie
fh1 = figure(1); clf
nr = 2 + double(iupdateR >= 1); % 3 rows for Pantuso model, 2 otherwise

% Large canvas so each of the nr-by-2 panels renders at a readable size.
% When saving (iSaveMov=1), the figure is also made invisible below, so a
% batch session renders the movie without opening a window. With iSaveMov=0
% the figure is shown on screen (a no-op under -nodisplay).
fh1.Position = [50 50 1800 1150];

Iunit = 1;    % Current in A
tunit = 1e6;  % Convert to us
LFS   = 24;
FS    = 20;
LW    = 2.4;

% Initialize file to save movie
if iSaveMov==1
    set(fh1,'Visible','off');
    % Initialize movie
    if ~exist('movName','var')
        movName = simCases{1};
    end
    % Saved under its own getTopResDir(movName) folder, same convention as
    % every other named output -- not necessarily the same folder as any
    % one input simulation, since this may overlay several of them.
    movResDir = getTopResDir(movName);
    if exist(movResDir,'dir')~=7
        mkdir(movResDir);
    end
    movPanels = VideoWriter([movResDir movName '.mp4'],'MPEG-4');
    movPanels.Quality = 95;
    open(movPanels);
    % getframe() on an invisible figure is unreliable (can return a blank
    % or stale bitmap depending on renderer) -- print() to a temp PNG and
    % read it back is the robust way to capture an off-screen figure.
    frameFile = [tempname '.png'];
end

% Main loop
for ksave=1:Ntsave

    for j=1:Nsims
        drs     = drsAll{j};
        Esave_j = EsaveAll{j};
        xE_j    = xEAll{j};

        subplot(nr,2,1) % Charge
        pQ = plot(drs.xV, drs.Qsave(ksave,:)*1e6, '-', 'LineWidth', LW);
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            hold off
            if Nsims>1
                legend(legName,'Location','NorthEast')
            end
            xlim([xmin xmax])
            ylim(Qlims*1e6)
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            xlabel('Distance (m)','FontSize',LFS)
            ylabel('Charge (\muC/m)','FontSize',LFS)
            title(['t = ' sprintf('%4.1f',tsave(ksave)*tunit) ' \mus'],'FontWeight','Normal');
            drawJunc();
        end

        subplot(nr,2,2) % Current
        pI = plot(drs.xI, drs.Isave(ksave,:)*Iunit, '-', 'LineWidth', LW);
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            hold off
            xlim([xmin xmax])
            ylim(Ilims*Iunit)
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            xlabel('Distance (m)','FontSize',LFS)
            ylabel('Current (A)','FontSize',LFS)
            [~, i10] = min(abs(drs.xI - 10));
            I10 = drs.Isave(ksave,i10);
            title(['I(x = 10 m) = ' sprintf('%1.1f',I10) ' A'],'FontWeight','Normal');
            drawJunc();
        end

        subplot(nr,2,3) % Electric field
        pE = semilogy(xE_j, max(abs(Esave_j(ksave,:)),Elims(1)*0.1), '-', 'LineWidth', LW);
        if j==1 && Nsims>1
            hold on
        end
        % The breakdown-field reference Ek (red dashed) is only meaningful as
        % "the threshold THIS run's front is crossing" -- with more than one
        % simulation overlaid it can't be attributed to either curve, so it is
        % only drawn for a single simulation.
        if iupdateR >= 1 && Nsims==1
            hold on
            semilogy(drs.xI, drs.Eksave(ksave,:), '--r', 'LineWidth', LW)
        end
        if j==Nsims
            hold off
            xlim([xmin xmax])
            ylim(Elims)
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on', ...
                'YTick',10.^(Elogmin:Elogmax))
            xlabel('Distance (m)','FontSize',LFS)
            ylabel('|E| (V/m)','FontSize',LFS)
            drawJunc();
        end

        subplot(nr,2,4) % Corona sheath radius rQ(x,t)
        prQ = plot(drs.xV, drs.rQsave(ksave,:), '-', 'LineWidth', LW);
        if j==1 && Nsims>1
            hold on
        end
        if j==Nsims
            % Gray dashed line: the electromagnetic signal horizon, x = xmax - v*t, i.e.
            % how far a signal launched from the tip at t = 0 can have travelled. It
            % leaves the plotted domain after (xmax-xmin)/v, about 14 us for a 3.3 km domain.
            line((xmax-tsave(ksave)*v)*[1 1],rQlims,'LineStyle','--','Color',.3*[1 1 1],'LineWidth',.9*LW)
            hold off
            xlim([xmin xmax])
            ylim(rQlims)
            set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
            xlabel('Distance (m)','FontSize',LFS)
            ylabel('r_Q (m)','FontSize',LFS)
            drawJunc();
        end

        if iupdateR >= 1
            Smin_j = drs.params.Smin;

            subplot(nr,2,5) % Temperature
            pT = plot(drs.xI, drs.Tsave(ksave,:)*1e-3, '-', 'LineWidth', LW);
            if j==1 && Nsims>1
                hold on
            end
            if j==Nsims
                hold off
                xlim([xmin xmax])
                ylim(Tlims*1e-3)
                set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on')
                xlabel('Distance (m)','FontSize',LFS)
                ylabel('Temperature (kK)','FontSize',LFS)
                drawJunc();
            end

            subplot(nr,2,6) % Conductivity
            pS = semilogy(drs.xI, max(drs.Ssave(ksave,:),Smin_j), '-', 'LineWidth', LW);
            if j==1 && Nsims>1
                hold on
            end
            if j==Nsims
                hold off
                xlim([xmin xmax])
                ylim(Slims)
                set(gca,'FontSize',FS,'TickDir','out','XMinorTick','on','YMinorTick','on', ...
                    'YTick',10.^(Slogmin:Slogmax))
                xlabel('Distance (m)','FontSize',LFS)
                ylabel('Conductivity (S/m)','FontSize',LFS)
                drawJunc();
            end
        end

    end
    drawnow;

    % Save movie?
    if iSaveMov==1
        print(fh1, '-dpng', '-r120', frameFile);
        writeVideo(movPanels, imread(frameFile));
    else
        pause(.0001)
    end

end

% Close file with movie
if iSaveMov==1
    close(movPanels);
    if exist(frameFile,'file'), delete(frameFile); end
end
