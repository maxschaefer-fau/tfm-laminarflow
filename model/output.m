%% output - evaluation of the output equation of the SSD to obtain concentration p 
% Calculates particle concentration p at point x = [r, phi, z] by the
% evaluation of the output equation of the SSD
%
% Syntax:  [out] = output(ybar, ftm, pickup, sim)
%
% Inputs:
%    ybar - system states (all orders)
%    ftm - stuct of tfm parameters 
%    pickup - pickup position
%    sim - struct of simulation parameters
%
% Outputs:
%    out - particle concentration at point observation point x
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

function [out] = output(ybar, ftm, pickup, sim)

p = zeros(ftm.N,length(sim.t));
cp = cell(1,ftm.M);

% exploit block partitioning - iterate over orders n 
for n = 1:ftm.N
    % calculate the output matrix for given order n
    cp{n} = outputMatrix_p(n, ftm.indexC{n}, ftm, pickup.pos);
    
    % distinguish between order n = 0 and n ~= 0 
    % exploits the relation between n > 0 and n < 0 (see
    % maximilianschaefer.org)
    if(ftm.n(n) ~= 0) 
        time_p = ybar{n};
        space_p = cp{n};
        p(n,:) = 2*real(space_p*time_p);
    else
        p(n,:) = cp{n}*ybar{n};
    end
end

out = sum(p,1);
end