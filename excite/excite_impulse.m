%% excite_impulse - initial values \bar{y}_init for an impulse 
%
% Syntax:  [excite, norm] = excite_impulse(ftm, index, n, excite_pos)
%
% Inputs:
%    ftm - struct containing all eigenvalues and eigenfunctions of tfm
%    index - index \mu for given order n 
%    n - current order n of Bessel functions
%    excite_pos - excitation position in the cylinder
%
% Outputs:
%    excite - initial values for the impulse
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

function [excite, norm] = excite_impulse(ftm, index, n, excite_pos)

r = excite_pos(1); 
phi = excite_pos(2); 
z = excite_pos(3); 

excite = zeros(size(index,1),1); 

for mu = 1:size(index,1)
    n_real = ftm.n(n); 
    
    m = index(mu,1);
    nu = index(mu,2); 
    
    l = ftm.lambda(nu); 
    k = ftm.k(n,m);
    
    excite(mu) = conj(ftm.kadj4(n_real, k, l, r, phi, z));
end

% normalize 
norm = 1;
end