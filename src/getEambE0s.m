function [E0, s] = getEambE0s(drs)
% getEambE0s  Ambient-field parameters of a saved DARTESS run, Eamb(x) = E0 + s*x
%
%   [E0, s] = getEambE0s(drs)
%
% E0 [V/m] is the ambient field at x = 0 (the conducting/cutoff boundary) and
% s [V/m^2] its gradient; s > 0 means the field weakens along the direction of
% dart leader propagation (-x). Reads every saved convention of the Eamb input:
%
%   dartess.m v0.102+     : drs.params.EambE0 / EambS (Eamb input = [E0 s])
%   dartess.m v0.97-v0.101: drs.params.Eamb = scalar, or [Val_at_xmax, Val_at_xmin]
%                           (a linear ramp between the run's own domain ends)
%   before v0.97          : no Eamb at all -> E0 = s = 0
%
% A 2-element drs.params.Eamb is ambiguous on its own, so the new form is
% recognized by the presence of drs.params.EambE0, never by the values.
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-10-09
% Last modification: 2026-10-09

p = drs.params;
if isfield(p, 'EambE0')
    E0 = p.EambE0;
    s  = p.EambS;
elseif isfield(p, 'Eamb') && ~isscalar(p.Eamb)
    xmx  = max(drs.xI);
    ldom = xmx - min(drs.xI);
    s    = (p.Eamb(1) - p.Eamb(2))/ldom;
    E0   = p.Eamb(1) - s*xmx;
elseif isfield(p, 'Eamb')
    E0 = p.Eamb;
    s  = 0;
else
    E0 = 0;
    s  = 0;
end
