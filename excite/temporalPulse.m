%% -------------- UNDER CONSTRUCTION -------------- %%

%% temporalPulse - temporal pulse shaping of external source fe
%
% Syntax:  fe_t = temporalPulse(excite_state_temp, diff, ftm, sim)
%
% Inputs:
%    excite_state_temp - flag to chose a temporal pulse shaping
%    diff - struct containing physical parameters
%    ftm - struct containing all eigenvalues and eigenfunctions of tfm
%    sim - struct containing simulation parameters
%
% Outputs:
%    fe_t - temporal pulse shaphing of external source fe
%       
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

function fe_t = temporalPulse(excite_state_temp, diff, ftm, sim)

switch excite_state_temp 
    
    case 'zero'
        fe_t = zeros(ftm.M*ftm.Nu, sim.len); 
        t0 = 0; 
        tr = 0;
    case 'dirac'
        t0 = 0.4; % position of dirac in [s]
        t0d = 0.1*Fs + 1;
        fe_t = zeros(ftm.M*ftm.Nu, sim.len);
        fe_t(:, t0d) = 1;
    
    case 'rect'    
        t0 = 0.1/diff.tau;   % pulse starting position [s]
        tr = 0.2/diff.tau;   % pulse width in [s]
        
        [fe_t] = excite_temp_rect(sim, t0, tr);
    case 'raised cosine' 
        t0 = 0.1/diff.tau;   % pulse starting position [s]
        tr = 0.3/diff.tau;   % pulse width in [s]
        
        [fe_t] = excite_temp_raised_cos(sim, t0, tr);
    otherwise
        warning('No temporal excitation function is set');
        return;
end

end