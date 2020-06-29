%% spatialPulse - initial conditions \bar{y}_init 
% Derives the expansion coefficients \bar{y}_init for the chosen spatial
% distribution of particles
%
% Syntax:  [yi, norm] = spatialPulse(initial_state, ftm, pos)
%
% Inputs:
%    initial_state - flag to chose an initial distribution
%    ftm - struct containing all eigenvalues and eigenfunctions of tfm
%    pos - spatial position of the initial distribution in the cylinder
%
% Outputs:
%    yi - initial conditions for the SSD
%    norm - normalization factor (used later to normalize output to the
%    number of injected particles)
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

function [yi, norm] = spatialPulse(initial_state, ftm, pos)

yi = cell(1,ftm.N);
norm = cell(1,ftm.N);

switch initial_state
    case 'zero'
        for n = 1:length(ftm.n)
           yi{n} = zeros(ftm.M*ftm.Nu,1);
           norm{n} = 0; 
        end
    case 'pointSingle'
        % pos: re, phie, ze, r0, phi0, z0
        for n = 1:length(ftm.n)
           [yi{n}, norm{n}] = ...
               excite_point(ftm, ftm.indexC{n}, n, pos); 
        end
        
    case 'circularSingle'
        % pos: re (doesnt matter), phi (doesnt matter), ze (center
        % position), r0, z0
        for n = 1:length(ftm.n)
           [yi{n}, norm{n}] = ...
               excite_circular(ftm, ftm.indexC{n}, n, pos); 
        end
        
    case 'impulse'
       % pos: re, phie, ze
       for n = 1:length(ftm.n)
          [yi{n}, norm{n}] = ...
              excite_impulse(ftm, ftm.indexC{n}, n, pos); 
       end
    otherwise
        warning('No initial conditions set.')
        return; 
end


end