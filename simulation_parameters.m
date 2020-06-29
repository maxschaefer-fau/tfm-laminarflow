%% simulation_parameters - Define basic simulation parameters
%
% Syntax:  [sim] = simulation_parameters(Fs, duration)
%
% Inputs:
%    Fs - sampling frequency
%    duration - normalized simulation duration
%
% Outputs:
%    sim - struct containing simulation parameters
%       sim.T - sampling time step
%       sim.len - length of simulation duration
%       sim.t - normalized simulation times 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
% Author: Dr.-Ing. Maximilian Schaefer, 
% University of Erlangen-Nuremberg
% email address: max.schaefer@fau.de
% Website: maximilianschaefer.org
% 28. June 2020; Last revision: 28. June 2020

function [sim] = simulation_parameters(Fs, duration)

sim.Fs = Fs;
sim.duration = duration;

sim.T = 1/sim.Fs;
sim.len = round(sim.Fs*sim.duration + 1);
sim.t = 0:sim.T:sim.duration;
end