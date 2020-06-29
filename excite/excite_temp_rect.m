%% -------------- UNDER CONSTRUCTION -------------- %%

%% excite_temp_rect - rect shaped 
% Creates a temporal pulse of length tr at starting position t0 shaped by a
% rectangular function. 
% Syntax:  [fe_t] = excite_temp_rect(sim, t0, tr)
%
% Inputs:
%    sim - struct containing simulation parameters
%    t0 - starting time of the pulse 
%    tr - width of the pulse
%
% Outputs:
%    fe_t - 1xsim.len array with temporal pulse of width tr at position t0
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
% March 2020; Last revision: 28. June 2020

function [fe_t] = excite_temp_rect(sim, t0, tr)
   
fe_t = zeros(1,sim.len);
rect = zeros(1,tr*sim.Fs);

for k = 1:length(rect)
    tf = sim.t(k);
    if(tf <= tr)
        rect(k) = 1;
    else
        rect(k) = 0;
    end
end

t1 = t0*sim.Fs + 1;
t2 = t1 + tr*sim.Fs -1;

fe_t(t1:t2) = rect;

end