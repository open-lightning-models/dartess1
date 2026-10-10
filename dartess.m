% DARTESS -- Dart-leader Acceleration & Retardation via Thermo-Electric
%           Simulation Scheme
%
% Geometry: horizontal dart leader in an intracloud flash
%
% DARTESS solves the Telegrapher's Equations along a partially-decayed
% intracloud channel, closed by a temperature- and field-dependent
% conductivity that is evolved in time, so the acceleration and retardation
% of the reactivation front are emergent rather than prescribed.
%
% Electrodynamics solver options (iQES flag):
%   iQES = 0 - full electrodynamics: telegrapher's equations (FDTD)
%   iQES = 1 - quasi-electrostatic: charge continuity + kernel-based
%              potential, I from Ohm's law
%
% Main state variables (capital letters):
%   I  - current [A]                     (at xI grid)
%   Q  - charge per unit length [C/m]    (at xV grid; primary state var)
%   R  - resistance per unit length [Ohm/m] (at xI grid)
%   C  - capacitance per unit length [F/m]  (at xV grid; per-node array derived from the
%                                         corona sheath radius rQ)
%   rQ - corona sheath radius [m]        (at xV grid; constant unless the sheath
%                                         expansion is enabled, see tauQss)
%   T  - temperature [K]                 (at xI grid, thermodynamic model only)
%   S  - electrical conductivity [S/m]   (at xI grid, thermodynamic model only)
%   V  - potential [V]                   (at xV grid; TOTAL potential: the fixed ambient
%                                         part Vamb plus the self-consistent part (Q-Q0)/C.
%                                         The uniform Q0 background is a static baseline)
%   E  - electric field [V/m]            (derived: E = -dV/dx, at xE grid; total field)
%   Vamb - ambient potential [V]         (at xV grid; fixed, built once from the optional
%                                         input Eamb; 0 if Eamb = 0, the default)
%
% Output: a structure drs saved to <getTopResDir(caseName)>/drs.mat, together
% with a copy of this script. See README.md for the list of saved fields.
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-06-05
% Last modified: 2026-10-09
% Code version: v1.0
%
clc; clear all

codeVersion = 'v1.0';   % saved to drs.version; keep in sync with the header above

% getTopResDir.m and the other helpers live in ./src, relative to the folder
% this script is run from.
if exist('src','dir'), addpath('src'); end

tic


%% Universal constants
e0 = 8.85e-12;
mu0 = 4*pi*1e-7;


%% Input params
% The values below are the reference case of the main article (Scenario I,
% Ilead = 100 A, sigma_coff = 1e-3 S/m).

% Spatial
xmin = -3000;       % Left (far) end of the domain [m]
xmax = +300;        % Right end of the domain = positive leader tip [m]
dx   = 0.5;         % Spatial step [m] (the QES disk kernel additionally assumes dx <= rQss)

% Temporal
tmax = 800e-6;      % Simulation duration [s] (the early stop at the domain edge may end the run sooner)
dt   = 10e-9;       % MAXIMUM time step [s]; the adaptive scheme steps below it as needed (at dx = 0.5 m the Courant limit, ~0.8 ns, binds)

% Trans. line params. (linear params)
% Channel junctions: entry 1 is the LEADER TIP and must sit at xmax; entry n injects the
% constant current Ijunc(n) at xjunc(n) and sets the initial conductivity of the section
% x <= xjunc(n), down to the next junction (the last one owns everything to xmin). A
% NEGATIVE S0junc(n) (e.g. -1) means that section is linearly ramped in x between its
% neighbors S0junc(n-1) and S0junc(n+1) instead of filled with a constant; entries 1 and
% Njunc can never be a ramp (no neighbor on one side).
% Uncomment exactly one of the three configurations below.

%--% Scenario I (conducting section followed by a cutoff section) -- reference case
xjunc  = [xmax,    0];  % Junction positions [m]; xjunc(1) must equal xmax (the leader tip)
Ijunc  = [ 1e2,    0];  % Injected current at each junction [A], +x sign (positive deposits negative charge); entry 1 is the leader tip current Ilead
S0junc = [ 1e2, 1e-3];  % Initial conductivity of the section starting at each junction [S/m]; a negative entry ramps linearly between its neighbors instead

% %--% Scenario II (step up in conductivity beyond x = -1000 m)
% xjunc  = [xmax,    0, -1000];
% Ijunc  = [ 1e2,    0,    0];
% S0junc = [ 1e2, 1e-3, 1e-2];

% %--% Scenario III (linear conductivity ramp from x = 0 to x = -1000 m)
% xjunc  = [xmax,    0, -0.5*dx, -1000];
% Ijunc  = [ 1e2,    0,       0,    0];
% S0junc = [ 1e2, 1e-3,      -1, 1e-1];

% Ambient electric field Eamb(x) = E0 + s*x, given as Eamb = [E0 s]: E0 [V/m] is the field at
% the cutoff point x = 0 and s [V/m^2] its gradient; s > 0 means the field weakens along the
% direction of dart leader propagation (-x). A scalar is a uniform field (s = 0). Positive
% values point toward xmax. The profile does not depend on xmin/xmax. 0 (default) = off.
Eamb = 0;    % [V/m] or [V/m, V/m^2]

Q0     = 0;           % Initial uniform charge per unit length [C/m]: residual charge on the decayed channel, superposed as a static (field-free) baseline
rI     = 3e-3;        % Current-carrying radius [m]: sets the ohmic cross section (pi*rI^2) and the inductance L
rQss   = 0.3;         % Steady-state corona sheath radius [m]: sets C (iQES=0) and the kernel disk radius (iQES=1)
lengthCL = 3e3;       % Reference channel length [m] used ONLY in the L/C log(length/radius) formulas -- independent of the domain (xmin/xmax)
tauQss = -1;          % Corona sheath expansion time constant [s]; 0 -> rQ = rI constant, <0 or >1e3*tmax -> rQ = rQss constant (no expansion)

% Electrodynamics solver
iQES = 0;           % =0 full electrodynamics (telegrapher's eqs, FDTD), =1 quasi-electrostatic (continuity + kernel potential)

% Adaptive time step (weights on the competing stability/accuracy limits)
wMDR  = 0.5;        % Weight on the Maxwell dielectric relaxation time tauMDR (iQES=1 electrodynamic limit)
wCFL  = 0.4;        % Weight on the Courant condition dx/v (iQES=0 electrodynamic limit)
wIon  = 0.1;        % Weight on 1/max(nu): fastest net Townsend ionization rate anywhere (iupdateR>=1 only)
wHeat = 0.1;        % Weight on min(T./max(dTdtHeat,dTdtCool)): fastest relative heating-or-cooling rate anywhere (iupdateR>=1 only)
dtmin = 1e-12;      % Adaptive-dt floor [s]: run stops with a warning if the combined limit drives dt below this

% Resistance model
iupdateR = 1;       % =0 R=const (frozen at its t=0 profile), >=1 thermodynamic (nonlinear) model

% Thermodynamic model params (iupdateR >= 1)
rg     = 5e-3;      % Thermal (gas) radius [m], rg>=rI
tauLTE = 10e-6;     % LTE relaxation time [s]
h      = 7;         % Altitude of dart leader above MSL [km]

% Townsend ionization/attachment params (Luque & Ebert 2009), sea-level values
nui0 = 1.65e4;      % Ionization rate prefactor mu0*alphai0 [m/(V s)]
nua0 = 76;          % Attachment rate prefactor mu0*alphaa0 [m/(V s)]
Ai   = 2e7;         % Ionization critical field [V/m]
Aa   = 3e6;         % Attachment critical field [V/m]
Ek0  = 3.16e6;      % Breakdown field [V/m] (matches (Ai-Aa)/log(nui0/nua0); set independently)

% Leader-current cutoff. Applies to the TIP current Ijunc(1) only; the other junctions
% inject constant current for the whole run.
Ioff = 1e-3;        % Floor current the tip current decays toward once the cutoff is active [A]
toff = 2*tmax;      % Time after which the tip current starts decaying [s]; toff >= tmax disables the cutoff

% Save (which time instants should be saved?)
Ntsave = 101; tsave = linspace(0,tmax,Ntsave); % Defined via number of instants to be saved
% tsave = [2 6 10]*1e-6; Ntsave = length(tsave);  % Alternative: defined via particular instants

% Flags
iSavVars  = 1; % =0 don't save vars, =1 does save it
caseName  = 'referenceCase'; % results folder name, under getTopResDir.m's root -- a run with the same caseName overwrites it


%% Secondary constants: TL electrical params
ldom = xmax - xmin;              % Total domain length [m] (sizes the solution-window buffer dibuf; NOT used for L/C, see lengthCL)
L  = mu0/2/pi*log(lengthCL/rI);  % Linear inductance [H/m] (constant: current flows within rI)

% Corona sheath radius (rQ) regime, resolved from tauQss. rQ itself is a per-xV-node
% array, set once the grid size is known (see "Independent vars"). Under QES the
% sheath radius is always the kernel radius rQss.
if iQES == 1
    iRQdynamic = 0;               % QES kernel is fixed at rQss
elseif tauQss == 0
    iRQdynamic = 0;                % no expansion: sheath stays at the current-carrying radius
elseif tauQss < 0 || tauQss > 1e3*tmax
    iRQdynamic = 0;                % already at steady state: constant radius rQss
else
    iRQdynamic = 1;                % per-node: each segment activates independently once triggered (see main loop)
end


%% Secondary constants: Plasma properties
delta_amb = exp(-h/9.7);          % Altitude air density reduction (Pantuso & da Silva 2024, Eq. 2)
Tamb      = 300;                  % Ambient temperature [K]
Smin      = 0;                    % Lower cap on conductivity [S/m] (numerical floor, not ambient air conductivity)
Tmax      = 40e3;                 % Cap = maximum possible, physically viable temperature


% Thermodynamic model secondary constants
if iupdateR >= 1
    % Thermodynamic coefficient functions (Pantuso & da Silva 2024)
    fC     = @(Tx) delta_amb*8.2e2*(Tx/300).^(-0.5);                    % Heat capacity [J/(m^3 K)]
    fkappa = @(Tx) 6.2e-2*(Tx/300);                                     % Thermal conductivity [W/(m K)]
    fSLTE  = @(Tx) 245*((Tx/300).^(1.125)).*exp(-13539./Tx);            % LTE conductivity [S/m]
    feps   = @(Tx) 4e11*((Tx/300).^(0.63)).*exp(-98000./Tx);            % Net emission coeff [W/m^3]
    fetaT  = @(Tx) 0.1 + 0.9*(tanh((Tx-1200)/300)+1)/2;                 % Joule heating efficiency, temperature factor [-]

    % T and S caps must be mutually consistent via S_LTE: Tmin solves S_LTE(Tmin) = Smin; Smax = S_LTE(Tmax).
    % Smin is capped first (at 1% of the lowest section conductivity), so Tmin matches the final Smin.
    Smin = min(Smin, 1e-2*min(S0junc(S0junc >= 0)));  % negative entries are ramp markers, not conductivities -- excluded here
    Tmin = fzero(@(x) fSLTE(x) - max(Smin, 7e-18), [Tamb Tmax]); Tmin = max(Tmin, Tamb);
    Smax = fSLTE(Tmax);
end


%% Check input and issue warnings
if ~(iupdateR == 0 || iupdateR >= 1)
    error('iupdateR must be 0 (constant R) or >=1 (thermodynamic model); got %g.', iupdateR);
end
% The Courant condition (iQES=0) and the Maxwell relaxation time (iQES=1) are
% enforced every step by the adaptive dt in the main loop.
if iQES == 1 && dx > rQss
    % Each segment is represented as a uniformly charged disk of radius rQss and
    % thickness dx (Pantuso & da Silva 2024, after their Eq. 6), which assumes dx <= rQss.
    % The TL closure C = 2*pi*e0/log(lengthCL/rQ) carries no such restriction.
    warning('QES disk kernel assumes dx <= rQss (dx = %g m, rQss = %g m).', dx, rQss);
end


%% Independent vars
xI   = xmax:-dx:xmin;           % decreasing: xI(1)=xmax (leader tip, right), xI(end)=xmin (left)
Nx   = length(xI);
[~, icoff] = min(abs(xI));      % first xI index in the originally-cutoff region (x <= 0)
ifront = icoff; %               % Starting index of  the wave front
dibuf  = round(5e-2*ldom/dx);   % Size of  buffer region ahead of dart leader wave (5% of the domain)
imaxq  = 1;                     % High-water mark of imax (QES kernel column window: charge beyond it never deviates from Q0)
xV   = xI(1:Nx-1) - dx/2;
xE   = xI(2:end-1);             % E-field grid: diff of V (at xV) maps to xI interior points

% Channel junctions. Each junction (i) injects the constant current Ijunc into the main
% channel -- a source in the continuity equation only, since the injected current arrives
% transversely and is a source of charge rather than an axial EMF -- and (ii) starts a new
% initial-conductivity section S0junc covering x <= xjunc, down to the next junction.
% Entry 1 is the leader tip, at xmax. Its current is applied through the boundary
% condition I(1) = -Ilead1 (equivalent to a qsrc entry at node 1, but it keeps I(1)/Igt
% meaningful as the channel-base current and lets Ioff/toff act on the tip alone), so it is
% excluded from qsrc below.
% Junctions are sorted rightmost-first, so "the next junction" means the next one leftward
% (xI/xV decrease with index); the tip (largest x) ends up first. Current is injected on
% the xV grid (where Q lives) while the conductivity step is applied on the xI grid (where
% S/T/R live), a half-cell offset inherent to the staggered grid.
Njunc = numel(xjunc);
if numel(Ijunc) ~= Njunc || numel(S0junc) ~= Njunc
    error('xjunc, Ijunc and S0junc must have the same number of elements (got %d, %d, %d).', ...
          Njunc, numel(Ijunc), numel(S0junc));
end
if Njunc < 1
    error('xjunc must have at least one entry: the leader tip at xmax.');
end
[xjunc, isrt]  = sort(xjunc, 'descend');   % tip first
Ijunc(1:Njunc)  = Ijunc(isrt);
S0junc(1:Njunc) = S0junc(isrt);
if xjunc(1) ~= xmax
    error(['xjunc(1) must be the leader tip at xmax = ' sprintf('%g', xmax) ...
           ' m (largest xjunc is ' sprintf('%g', xjunc(1)) ' m).']);
end
if any(diff(xjunc) == 0)
    error('xjunc entries must be distinct (two junctions at the same position).');
end
if min(xjunc) < xmin
    error(['xjunc entries must lie within the domain [' sprintf('%g, %g', xmin, xmax) '] m.']);
end

% Ramp sections: a negative S0junc(n) means the section is linearly ramped in x from the
% conductivity of the section to its right (S0junc(n-1)) to that of the section to its
% left (S0junc(n+1)), continuous with both. Both neighbors must be real (non-negative)
% conductivities.
for n = 1:Njunc
    if S0junc(n) < 0
        if n == 1 || n == Njunc
            error(['S0junc(' num2str(n) ') is negative (a ramp section), but ramp sections need a real ' ...
                   'conductivity on both sides -- entry 1 (the tip section) and entry ' num2str(Njunc) ...
                   ' (the last section) can never be a ramp.']);
        end
        if S0junc(n-1) < 0 || S0junc(n+1) < 0
            error(['S0junc(' num2str(n) ') is a ramp section, but its immediate neighbor is too ' ...
                   '(entry ' num2str(n-1) ' or ' num2str(n+1) ') -- chained ramp sections are not supported.']);
        end
    end
end

% Ambient electric field input validation
if ~(isnumeric(Eamb) && isvector(Eamb) && (numel(Eamb) == 1 || numel(Eamb) == 2))
    error('Eamb must be a scalar or a 2-element [E0 s] array (got %d elements).', numel(Eamb));
end
EambE0 = Eamb(1);                       % ambient field at x = 0 [V/m]
EambS  = 0; if numel(Eamb) == 2, EambS = Eamb(2); end   % ambient field gradient [V/m^2]

jjunc = zeros(1,Njunc);          % xV-grid index of each junction
qsrc  = zeros(1,Nx-1);           % lateral charge source [A/m] on the xV grid: dQ/dt gets -qsrc
for n = 1:Njunc
    [~, jjunc(n)] = min(abs(xV - xjunc(n)));       % nearest xV node (Q's grid)
    if n > 1                                        % entry 1 is the BC, NOT a qsrc entry
        qsrc(jjunc(n)) = qsrc(jjunc(n)) + Ijunc(n)/dx;  % accumulate (harmless if two ever share a cell)
    end
end
disp(['Junctions at x = ' sprintf('%g ', xV(jjunc)) 'm (snapped from ' sprintf('%g ', xjunc) 'm); tip current = ' sprintf('%g', Ijunc(1)) ' A']);

% Leftmost injecting junction (entries 2..N with a nonzero current). These are live from
% t = 0, so the active domain must reach past them (see the imax update in the main loop).
% Junctions that only mark a section boundary (Ijunc = 0) inject nothing and are not
% included. Set to 0 when nothing injects.
isrcj = false(1,Njunc); isrcj(2:end) = (Ijunc(2:end) ~= 0);
if any(isrcj)
    jsrcmax = max(jjunc(isrcj));
else
    jsrcmax = 0;
end

% Corona sheath radius (rQ): per-xV-node array. In the static regimes every node keeps
% the same constant value; in the dynamic regime every node starts at rI and grows once
% it locally crosses the breakdown-field trigger (see the end of the main loop). Under
% QES the sheath radius always equals the kernel radius rQss.
if iQES == 1
    rQ = rQss*ones(1,Nx-1);            % QES: pinned to the kernel radius
elseif iRQdynamic
    rQ = rI*ones(1,Nx-1);              % starts at rI everywhere; individual nodes grow once triggered
    isRQactive = false(1,Nx-1);        % one-way latch: true once a node has crossed the trigger, stays true
elseif tauQss == 0
    rQ = rI*ones(1,Nx-1);              % no expansion: frozen at rI
else
    rQ = rQss*ones(1,Nx-1);            % tauQss<0 or >1e3*tmax: frozen at rQss
end
C = 2*pi*e0./log(lengthCL./rQ);  % Linear capacitance [F/m], per xV node, from the sheath radius rQ
v = 1/sqrt(L*C(1));               % Speed of EM waves guided by the TL [m/s], evaluated at rQ(t=0)
Z = sqrt(L/C(1));                 % Impedance [Ohm], evaluated at rQ(t=0)

% Time array, preallocated for the expected number of steps: the input dt is only a ceiling,
% and on the TL branch the Courant limit usually sets the step (documentation.md, Sec. 6).
% The per-step arrays still grow if dt adapts below this estimate.
if iQES == 0
    dtEst = min(dt, wCFL*dx/v);
else
    dtEst = dt;
end
Nt = ceil(tmax/dtEst) + 1;
t  = zeros(1,Nt);                 % t(1) = 0; t(k+1) = t(k) + dt is set step by step in the main loop

% Ambient potential Vamb: closed-form antiderivative of the ambient field
% Eamb(x) = E0 + s*x, on the xV grid. With u = xmax - x (so E = -dV/dx = dV/du),
% Eamb = Etip - s*u, Etip = E0 + s*xmax, and Vamb = Etip*u - s*u^2/2.
% Reference Vamb = 0 at the tip (x = xmax, u = 0); only differences of Vamb enter
% the physics. Built once: a fixed background, not a state variable. Vamb is folded into
% V wherever V is rebuilt from Q, so V is the total potential, E = -dV/dx the total field,
% and the saved Vsave the total potential. At Eamb = 0, Vamb is identically 0.
u = xmax - xV;                          % distance from the tip, along xV; >= 0
Vamb = (EambE0 + EambS*xmax)*u - EambS*u.^2/2;   % s = 0: linear ramp (uniform field)


%% Preallocation
% Model vars
Q          = Q0*ones(1,Nx-1);   % Charge per unit length [C/m]; primary (TOTAL) state var, initialized uniform at Q0
V          = (Q-Q0)./C + Vamb;  % TOTAL potential (= Vamb at t=0, since Q = Q0)
I          = zeros(1,Nx);

% Initial conductivity profile (xI grid), from which R -- and, under the thermodynamic
% model, S and T -- are derived. Junction n owns every node with x <= xjunc(n) down to the
% next junction leftward; the last one owns everything down to xmin. A node exactly AT a
% junction belongs to the section that starts there, so with a junction at x = 0 the node
% at x = 0 is in the cutoff section. Ramp sections get a linear-in-x fill, continuous
% with S0junc(n-1) at their right edge and S0junc(n+1) at their left edge.
S0vec = zeros(1,Nx);
for n = 1:Njunc
    if n < Njunc
        msec = (xI <= xjunc(n)) & (xI > xjunc(n+1));
    else
        msec = (xI <= xjunc(n));
    end
    if S0junc(n) < 0
        slope = (S0junc(n+1) - S0junc(n-1)) / (xjunc(n+1) - xjunc(n));
        S0vec(msec) = S0junc(n-1) + slope*(xI(msec) - xjunc(n));
    else
        S0vec(msec) = S0junc(n);
    end
end
if any(S0vec == 0)
    error('%d grid nodes were not covered by any junction section -- check xjunc.', sum(S0vec==0));
end
R = 1 ./ (S0vec * pi * rI^2);   % Initial resistance per unit length [Ohm/m]
Qo         = zeros(1,Nx-1);
Io         = zeros(1,Nx);
Ro         = zeros(1,Nx);
% Save vars
Vsave    = zeros(Ntsave,Nx-1);
Qsave    = zeros(Ntsave,Nx-1);
Isave    = zeros(Ntsave,Nx);
Rsave    = zeros(Ntsave,Nx);
rQsave   = zeros(Ntsave,Nx-1);  % corona sheath radius [m] per xV node, at each snapshot (uniform unless the sheath expands)
Igt      = zeros(1,Nt);
xIwavpos = zeros(1,Nt);  % wave front position [m]; starts at x=0, holds previous value when undetected
Ewavpos  = zeros(1,Nt);  % |E| at the wave front position [V/m]; held in sync with xIwavpos

% Thermodynamic model arrays (iupdateR >= 1)
if iupdateR >= 1
    % Initial T consistent with the initial conductivity: solve S_LTE(T) = S0 once per
    % distinct conductivity in S0vec, then map back onto the grid. The brackets are the
    % full [Tamb Tmax] range (S_LTE is monotonic in T); the max(.,7e-18) floor keeps fzero
    % well-posed for vanishingly small S0. Conductivities above Smax are out of range.
    if max(S0vec) > Smax
        error(['Initial conductivities must not exceed Smax = S_LTE(Tmax) = ' ...
               sprintf('%.4g', Smax) ' S/m (got max(S0junc) = %.4g).'], max(S0vec));
    end
    S0uniq = unique(S0vec);
    T      = zeros(1, Nx);
    for n = 1:numel(S0uniq)
        T(S0vec == S0uniq(n)) = fzero(@(x) fSLTE(x) - max(S0uniq(n), 7e-18), [Tamb Tmax]);
    end
    S = S0vec;

    % LTE temperature of each section, saved to drs.params for reference. Ramp sections
    % have no single representative conductivity, so they get NaN.
    Tjunc = zeros(1,Njunc);
    for n = 1:Njunc
        if S0junc(n) < 0
            Tjunc(n) = NaN;
        else
            Tjunc(n) = fzero(@(x) fSLTE(x) - max(S0junc(n), 7e-18), [Tamb Tmax]);
        end
    end

    % Override R from S
    R = 1 ./ (S .* pi .* rI^2);

    % Save arrays for thermodynamic diagnostics
    Tsave  = zeros(Ntsave, Nx);
    Ssave  = zeros(Ntsave, Nx);
    Eksave = zeros(Ntsave, Nx);  % Effective breakdown field Ek0*delta [V/m] on xI grid

    % End-of-previous-step rates (all zero at t = 0, where E = 0): read by the adaptive
    % dt at the top of each iteration, refreshed by the corrector each step
    nu      = zeros(1, Nx);   % net Townsend ionization rate [1/s]
    dTdtHeat = zeros(1, Nx);   % Joule heating rate [K/s]
    dTdtCool = zeros(1, Nx);   % conduction + radiative cooling rate [K/s]
end

%% QES kernel: inverse-capacitance matrix (iQES = 1)
% Each xV segment is modeled as a uniformly charged disk of radius rQss and thickness dx
% (on-axis potential of a charged disk), in free space -- no image charges, since the
% horizontal intracloud channel has no ground plane. Built once at startup.
if iQES == 1
    kC = 1/(4*pi*e0);                                 % Coulomb's constant [N m^2/C^2]
    [Xj, Xi] = meshgrid(xV, xV);                      % Xi: observation segments (rows), Xj: source segments (columns)
    rk = abs(Xi - Xj);                                % axial distance between segment centers [m]
    Aker = (2*kC/rQss^2)*(sqrt(rk.^2 + rQss^2) - rk); % potential at Xi per unit charge at Xj [1/F]
    clear Xj Xi rk
    % The uniform Q0 background is a static baseline left out of V (as in the TL branch):
    % V is the potential of the deviation Q - Q0 alone, plus Vamb.
    % Electrodynamic stability (Maxwell dielectric relaxation time, Pantuso & da Silva
    % 2024, Eq. 19) is enforced every step by the adaptive dt in the main loop.
end


%% Main loop
% Snapshot 1 is the initial condition, saved before the first step (when tsave starts at
% t(1), the linspace default). The in-loop save runs at the end of each step.
ksave = 1;
saveIntFlag = 1;
if tsave(1) <= t(1)
    Isave(ksave,:) = -I;
    Vsave(ksave,:) = V;
    Qsave(ksave,:) = Q;
    Rsave(ksave,:) = R;
    rQsave(ksave,:) = rQ;
    if iupdateR >= 1
        Tsave(ksave,:)  = T;
        Ssave(ksave,:)  = S;
        Eksave(ksave,:) = Ek0 * delta_amb * Tamb ./ T;
    end
    tsave(ksave) = t(1);
    ksave = ksave + 1;
    if ksave > Ntsave
        saveIntFlag = 0;
        ksave = Ntsave;
    end
end
Ilead1 = Ijunc(1);    % tip current (junction entry 1), applied via the boundary condition
dt0 = dt;       % input (maximum) time step; dt adapts below it every step
k = 1;
while t(k) <= tmax

    % Global adaptive time step: dt is the smallest of dt0 and the active physical/
    % numerical time scales, each scaled by its own safety weight. Electrodynamics
    % contributes one term (set by iQES); the two thermodynamic terms are active when the
    % thermodynamic model runs. All candidates use the previous step's state.
    if iQES == 1
        tauMDR = pi*e0*rQss^2*max(min(R), eps); % Maxwell dielectric relaxation time (Pantuso & da Silva 2024, Eq. 19), via R = 1/(S*pi*rI^2)
        dt = min(dt0, wMDR*tauMDR);             % eps floor guards against dt collapsing to 0
    else
        dt = min(dt0, wCFL*dx/v);          % Courant condition (one-cell EM wave-crossing time)
    end
    if iupdateR >= 1
        dt = min(dt, wIon *(1/max(max(nu), eps)));                    % fastest net ionization rate anywhere (eps floor: guards against nu = 0 everywhere)
        dt = min(dt, wHeat*min(T./max(max(dTdtHeat,dTdtCool), eps))); % fastest relative heating-or-cooling rate anywhere (eps floor: same)
    end

    if dt < dtmin
        warning('Adaptive dt fell below dtmin (dt = %.2g s at t = %.3g s). Stopping.', dt, t(k));
        k = max(k-1, 1);   % step k has no state/diagnostics yet; trim to the last completed step
        break
    end
    if mod(k,1e5) == 0
        disp(['t = ' sprintf('%.4g',t(k)*1e6) ' us, dt = ' sprintf('%.3g',dt*1e9) ' ns, k = ' num2str(k)]);
    end

    % Active solution window: nodes 1..imax are updated this step
    if iQES == 0
        imax = min(ceil(t(k)*v/dx),Nx-1);  % TL: turn-on front travels at the EM wave speed v
        imax = min(imax, ifront+dibuf);
    else
        imax = min(ifront+dibuf, Nx-1);    % QES: potential responds instantaneously; window tracks the physical front
    end
    % Injecting junctions and a nonzero ambient field act from t = 0, so the window must
    % reach past them from the first step. Beyond the leftmost injecting junction it grows
    % at the EM signal speed v (the causal bound on how far the source's influence can have
    % reached); a fixed buffer would act as a reflecting wall, since the flux at cell imax
    % uses I(imax+1), which stays 0. Under QES there is no propagation delay, so the whole
    % domain is opened at once.
    if jsrcmax > 0 || any(Eamb ~= 0)
        if iQES == 0
            imax = max(imax, jsrcmax + max(dibuf, ceil(t(k)*v/dx)));
        else
            imax = Nx-1;
        end
        imax = min(imax, Nx-1);
    end
    if iQES == 1
        imaxq = max(imaxq, imax);          % high-water mark (kernel column window; ifront can transiently regress)
    end

    if t(k)>toff
       Ilead1 = Ioff + (Ijunc(1) - Ioff)*exp(-(t(k)-toff)/10e-6);
    end

    % Store info
    Io = I;
    Qo = Q;
    Ro = R;

    % Save thermodynamic state at start of step (needed for trapezoidal corrector)
    if iupdateR >= 1
        To = T; So = S;
    end

    % First half step
    if iQES == 0
        aR = dt*R/2/L;
        for i=2:imax
            I(i) = I(i)*(1-aR(i))/(1+aR(i)) - dt/L*(V(i)-V(i-1))/dx/(1+aR(i));
            Q(i) = Q(i) - dt*(I(i+1)-I(i))/dx - dt*qsrc(i);
            V(i) = (Q(i)-Q0)/C(i) + Vamb(i);
        end
        I(1) = -Ilead1;
        Q(1) = Q(1) - dt*(I(2)-I(1))/dx - dt*qsrc(1);
        V(1) = (Q(1)-Q0)/C(1) + Vamb(1);
    else
        jj = 2:imax;
        I(jj) = -(V(jj)-V(jj-1))/dx ./ R(jj);   % Ohm's law I = E/R (array-native sign convention)
        I(1) = -Ilead1;
        iq = 1:imax;
        Q(iq) = Q(iq) - dt*(I(iq+1)-I(iq))/dx - dt*qsrc(iq);
        V = (Aker(:,1:imaxq)*(Q(1:imaxq)-Q0)')'*dx + Vamb;   % kernel potential of the deviation (charge beyond imaxq stays at Q0) plus the ambient background
    end
    % Thermodynamic predictor - rates at (To, So), advance by dt/2
    if iupdateR >= 1
        E = zeros(1, Nx);
        if imax >= 2
            jj = 2:min(imax, Nx-1);
            E(jj) = -(V(jj) - V(jj-1)) / dx;   % total field (V is total)
        end
        E(1) = E(2); E(Nx) = E(Nx-1);

        % Equation coefficients, f = f(T or E)
        Cheat  = delta_amb*8.2e2*(To/300).^(-0.5);
        kappa  = 6.2e-2*(To/300);
        SLTE   = 245*((To/300).^(1.125)).*exp(-13539./To);
        epsNEC = 4e11*((To/300).^(0.63)).*exp(-98000./To);
        eta    = 0.1 + 0.9*(tanh((To-1200)/300)+1)/2;
        delta  = delta_amb * Tamb ./ To;
        nui    = nui0 * abs(E) .* exp(-Ai*delta./abs(E));
        nua    = nua0 * abs(E) .* exp(-Aa*delta./abs(E));
        nu     = max(nui - nua, 0);

        % Update equations' rhs
        dTdtHeat = eta.*So.*E.^2 ./ Cheat;                          % Joule heating rate [K/s]
        dTdtCool = ((4*kappa/rg^2).*(To-Tamb) + 4*pi*epsNEC) ./ Cheat; % conduction + radiative cooling rate [K/s]
        dTdt_1    = dTdtHeat - dTdtCool;
        dSdt_1    = nu.*So - (So-SLTE)/tauLTE;

        T     = min(max(To + (dt/2)*dTdt_1, Tmin), Tmax);
        S = min(max(So + (dt/2)*dSdt_1, Smin), Smax);

        R = 1 ./ (S .* pi .* rI^2);
    end

    % Second half step
    if iQES == 0
        aR = dt*(R + Ro)/4/L;
        for i=2:imax
            I(i) = Io(i)*(1-aR(i))/(1+aR(i)) - dt/L*(V(i)-V(i-1))/dx/(1+aR(i));
            Q(i) = Qo(i) - dt*(I(i+1)-I(i))/dx - dt*qsrc(i);
            V(i) = (Q(i)-Q0)/C(i) + Vamb(i);
        end
        I(1) = -Ilead1;
        Q(1) = Qo(1) - dt*(I(2)-I(1))/dx - dt*qsrc(1);
        V(1) = (Q(1)-Q0)/C(1) + Vamb(1);
    else
        jj = 2:imax;
        I(jj) = -(V(jj)-V(jj-1))/dx*2./(R(jj)+Ro(jj));   % Ohm's law with half-step-averaged R (mirrors the TL's aR)
        I(1) = -Ilead1;
        iq = 1:imax;
        Q(iq) = Qo(iq) - dt*(I(iq+1)-I(iq))/dx - dt*qsrc(iq);
        V = (Aker(:,1:imaxq)*(Q(1:imaxq)-Q0)')'*dx + Vamb;
    end
    % Thermodynamic corrector - rates at (T_half, S_half), trapezoidal update from To
    if iupdateR >= 1
        E = zeros(1, Nx);
        if imax >= 2
            jj = 2:min(imax, Nx-1);
            E(jj) = -(V(jj) - V(jj-1)) / dx;   % total field (V is total)
        end
        E(1) = E(2); E(Nx) = E(Nx-1);

        % Equation coefficients, f = f(T or E)
        Cheat  = delta_amb*8.2e2*(T/300).^(-0.5);
        kappa  = 6.2e-2*(T/300);
        SLTE   = 245*((T/300).^(1.125)).*exp(-13539./T);
        epsNEC = 4e11*((T/300).^(0.63)).*exp(-98000./T);
        eta    = 0.1 + 0.9*(tanh((T-1200)/300)+1)/2;
        delta  = delta_amb * Tamb ./ T;
        nui    = nui0 * abs(E) .* exp(-Ai*delta./abs(E));
        nua    = nua0 * abs(E) .* exp(-Aa*delta./abs(E));
        nu     = max(nui - nua, 0);

        % Update equations' rhs
        dTdtHeat = eta.*S.*E.^2 ./ Cheat;                          % Joule heating rate [K/s]; also used by the next step's adaptive dt
        dTdtCool = ((4*kappa/rg^2).*(T-Tamb) + 4*pi*epsNEC) ./ Cheat; % conduction + radiative cooling rate [K/s]
        dTdt_2    = dTdtHeat - dTdtCool;
        dSdt_2    = nu.*S - (S-SLTE)/tauLTE;

        % Trapezoidal update from To, So
        T = min(max(To + (dt/2)*(dTdt_1+dTdt_2), Tmin), Tmax);
        S = min(max(So + (dt/2)*(dSdt_1+dSdt_2), Smin), Smax);

        R = 1 ./ (S .* pi .* rI^2);
    end

    % Wave front: location of max |E| within the originally-cutoff region (x <= 0),
    % searched only over a window that follows the front (out to ifront+dibuf). This
    % keeps the tracker on the propagating front rather than on distant |E| peaks, such as
    % one raised from t = 0 by a current-injecting junction ahead of the front. Until the
    % active domain reaches x = 0, the previous position is held (xIwavpos starts at 0).
    % The front position is not allowed to move backward (toward +x).
    iEnd = min(imax, ifront + dibuf);
    if iEnd >= icoff
        Eact = abs(-(V(icoff:iEnd) - V(icoff-1:iEnd-1)) / dx);  % |E| at xI(icoff:iEnd); total field
        [Emax, imaxE] = max(Eact);
        ifront = icoff - 1 + imaxE;
        xIwavpos(k) = xI(ifront);
        Ewavpos(k)  = Emax;
    elseif k > 1
        xIwavpos(k) = xIwavpos(k-1);
        Ewavpos(k)  = Ewavpos(k-1);
    end
    if k>1 && xIwavpos(k) > xIwavpos(k-1)
        xIwavpos(k) = xIwavpos(k-1);
        Ewavpos(k)  = Ewavpos(k-1);
    end

    % Save variables as a function of both space and time
    if t(k)>=tsave(ksave) && saveIntFlag==1

        Isave(ksave,:) = -I;
        Vsave(ksave,:) = V;
        Qsave(ksave,:) = Q;
        Rsave(ksave,:) = R;
        rQsave(ksave,:) = rQ;
        if iupdateR >= 1
            Tsave(ksave,:)   = T;
            Ssave(ksave,:)   = S;
            Eksave(ksave,:) = Ek0 * delta_amb * Tamb ./ T;
        end


        % Update ksave
        ksave = ksave + 1;
        if ksave>Ntsave
            saveIntFlag = 0;
            ksave = Ntsave;
        end
    end

    % Save channel-base current as a function of time
    Igt(k) = I(1);

    % Corona sheath radius (rQ) dynamics (iRQdynamic only): a node is activated once the
    % radial field at its sheath boundary, Er_Q (Gauss's law, from the deviation charge
    % Q-Q0), exceeds the local breakdown field Ek_Q; active nodes then relax toward rQss
    % semi-implicitly (Crank-Nicolson), once per full step. Ek_Q is Ek0*delta averaged
    % from the xI grid (where T lives) onto xV (where Q/rQ live), or an altitude-only value
    % if the thermodynamic model is off.
    if iRQdynamic
        if iupdateR >= 1
            delta_I = delta_amb * Tamb ./ T;
            Ek_Q = Ek0 * (delta_I(1:Nx-1) + delta_I(2:Nx)) / 2;   % xI -> xV
        else
            Ek_Q = Ek0 * delta_amb * ones(1,Nx-1);                      % no T model: altitude-only value
        end
        Er_Q = abs(Q-Q0) ./ (2*pi*e0*rQ);
        isRQactive = isRQactive | (Er_Q > Ek_Q);   % one-way latch: once triggered, always active

        bQ = dt/(2*tauQss);
        rQnew = ((1-bQ)*rQ + 2*bQ*rQss)/(1+bQ);
        rQ(isRQactive) = rQnew(isRQactive);          % inactive nodes stay exactly at rI
        C = 2*pi*e0./log(lengthCL./rQ);
    end

    % Stop early once the wave front is within dibuf cells of the left edge (the solution
    % window would otherwise be clipped by xmin)
    if (Nx - ifront) <= dibuf
        break
    end

    % Advance the clock (t grows past its initial allocation if dt shrank)
    t(k+1) = t(k) + dt;
    k = k + 1;

end
if t(k) > tmax
    k = k - 1;   % normal exit: step k was never executed (breaks leave t(k) <= tmax, so they are unaffected)
end


%% Trim time series to the steps actually taken (with adaptive dt the final step
%% count is not known a priori, so a final snapshot save is always made)
if k > 0

    % Force a final save at the instant the loop was stopped
    Isave(ksave,:) = -I;
    Vsave(ksave,:) = V;
    Qsave(ksave,:) = Q;
    Rsave(ksave,:) = R;
    rQsave(ksave,:) = rQ;
    if iupdateR >= 1
        Tsave(ksave,:)  = T;
        Ssave(ksave,:)  = S;
        Eksave(ksave,:) = Ek0 * delta_amb * Tamb ./ T;
    end
    tsave(ksave) = t(k);
    Ntsave = ksave;

    % Trim time-stepped vars to the instant the loop was stopped
    t        = t(1:k);
    Igt      = Igt(1:k);
    xIwavpos = xIwavpos(1:k);
    Ewavpos  = Ewavpos(1:k);

    % Trim snapshot vars to the instants actually saved
    tsave = tsave(1:Ntsave);
    Vsave = Vsave(1:Ntsave,:);
    Qsave = Qsave(1:Ntsave,:);
    Isave = Isave(1:Ntsave,:);
    Rsave = Rsave(1:Ntsave,:);
    rQsave = rQsave(1:Ntsave,:);
    if iupdateR >= 1
        Tsave   = Tsave(1:Ntsave,:);
        Ssave   = Ssave(1:Ntsave,:);
        Eksave = Eksave(1:Ntsave,:);
    end
end


%% Variable ranges (saved to drs; gen1Dmov.m reads drs.Ilims)
Rmin = min(Rsave(:));
if Rmin==0
    Rmin = -1;
else
    Rmin = 0;
end
Rmax = max(Rsave(:));
if Rmax==0
    Rmax = 1;
else
    Rmax = 1.2*Rmax;
end
Rlims = [Rmin Rmax];

Vmin = min(Vsave(:));
Vmax = max(Vsave(:));
dV = 0.05*(Vmax-Vmin);
Vlims = [Vmin-dV Vmax+dV];

Imin = min(Isave(:));
Imax = max(Isave(:));
dI = 0.05*(Imax-Imin);
Ilims = [Imin-dI Imax+dI];


%% Save variables (drs is structure with important vars)
if iSavVars>0

    % Each run gets its own folder, named after caseName, under the top-level results
    % directory (getTopResDir.m). An existing folder with the same name is cleared first.
    topResDir = getTopResDir(caseName);
    if exist(topResDir,'dir')==7
        delete([topResDir '*']);
    else
        mkdir(topResDir);
    end

    % Code version
    drs.version = codeVersion;

    % Save input parameters
    drs.params.xmin     = xmin;
    drs.params.xmax     = xmax;
    drs.params.tmax     = tmax;
    drs.params.dt       = dt0;   % input (maximum) time step; per-step dt history = diff(drs.t)
    drs.params.dx       = dx;
    drs.params.rI       = rI;
    drs.params.S0junc   = S0junc;  % initial conductivity [S/m] of the section starting at each junction (negative = ramp section); the realized profile is drs.Ssave(1,:)
    drs.params.L        = L;
    drs.params.C        = C;      % final value (C(rQ) evolves if the sheath expands -- see drs.rQsave for the full history)
    drs.params.Q0       = Q0;
    drs.params.Eamb     = Eamb;   % ambient field input as given: scalar [V/m] or [E0 s] [V/m, V/m^2] (v0.102+; v0.97-v0.101 used [Val_at_xmax, Val_at_xmin]); 0 = off
    drs.params.EambE0   = EambE0; % ambient field at x = 0 [V/m], Eamb(x) = EambE0 + EambS*x
    drs.params.EambS    = EambS;  % ambient field gradient [V/m^2]
    drs.params.v        = v;      % evaluated at rQ(t=0)
    drs.params.Z        = Z;      % evaluated at rQ(t=0)
    drs.params.rQss     = rQss;
    drs.params.lengthCL = lengthCL;  % reference length [m] used in the L/C log(length/radius) formulas
    drs.params.tauQss   = tauQss;
    drs.params.xjunc    = xjunc;   % junction positions [m], sorted rightmost-first; entry 1 is the leader tip at xmax
    drs.params.Ijunc    = Ijunc;   % injected current at each junction [A]; entry 1 is the tip current (applied as the BC)
    drs.params.jjunc    = jjunc;   % xV-grid index of each junction
    drs.params.xjuncSnapped = xV(jjunc);   % positions actually used (nearest xV node)
    drs.params.Ioff     = Ioff;
    drs.params.toff     = toff;
    drs.params.wMDR     = wMDR;
    drs.params.wCFL     = wCFL;
    drs.params.wIon     = wIon;
    drs.params.wHeat    = wHeat;
    drs.params.dtmin    = dtmin;

    if iupdateR >= 1
        drs.params.rg        = rg;
        drs.params.tauLTE   = tauLTE;
        drs.params.Smin = Smin;
        drs.params.Smax = Smax;
        drs.params.Tmin = Tmin;
        drs.params.Tmax = Tmax;
        drs.params.Tamb      = Tamb;

        drs.params.h         = h;
        drs.params.delta_amb = delta_amb;
        drs.params.Ek0       = Ek0;
        drs.params.Tjunc     = Tjunc;  % LTE temperature of each junction's section [K]

        drs.params.nui0      = nui0;
        drs.params.nua0      = nua0;
        drs.params.Ai        = Ai;
        drs.params.Aa        = Aa;
    end

    % Save flags
    drs.flags.iQES      = iQES;
    drs.flags.iupdateR  = iupdateR;
    drs.flags.iSavVars  = iSavVars;
    drs.flags.caseName  = caseName;

    % Save variables
    drs.t = t;
    drs.tsave = tsave;
    drs.xI = xI;
    drs.xV = xV;
    drs.Vamb = Vamb;   % ambient potential on the xV grid [V]; Vsave is the TOTAL potential, so the self-consistent part alone is Vsave - Vamb (per row)
    drs.Vsave = Vsave;
    drs.Qsave = Qsave;
    drs.Isave = Isave;
    drs.Rsave = Rsave;
    drs.rQsave = rQsave;
    drs.Igt = Igt;
    if iupdateR >= 1
        drs.Tsave     = Tsave;
        drs.Ssave = Ssave;
        drs.Eksave    = Eksave;
        drs.T     = T;
        drs.S = S;
    end
    drs.xIwavpos = xIwavpos;
    drs.Ewavpos  = Ewavpos;

    % Save var limits
    drs.Vlims = Vlims;
    drs.Ilims = Ilims;
    drs.Rlims = Rlims;

    save([topResDir 'drs'],'drs');

    % Save a copy of this script exactly as it was at run time (whatever its name)
    copyfile([mfilename('fullpath') '.m'], topResDir);
end

toc
