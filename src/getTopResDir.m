function fullpath = getTopResDir(caseName)
%
% fullpath = getTopResDir(caseName)
%
% Returns the full path to a DARTESS simulation's results folder,
% <top-level results directory>/<caseName>/. Every run is distinguished only
% by its caseName, so a run with an existing caseName overwrites it.
%
% With no argument, returns the top-level directory itself.
%
% The default top-level directory is ./runs/ inside the folder that holds
% dartess.m (one level above this file). Change topresdir below to keep
% results elsewhere.
%
% Author: Caitano L. da Silva, NMT
% Created: 2026-08-09
%

%% Top results directory (change as needed)
topresdir = [fileparts(fileparts(mfilename('fullpath'))) filesep 'runs' filesep];

%% Output path
if nargin==0 || isempty(caseName)
    fullpath = topresdir;
else
    fullpath = [topresdir caseName '/'];
end

end
