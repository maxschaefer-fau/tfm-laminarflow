%% outputMatrix_p - calculates output matrix c_1 for order n
% calculates the output matrix c_1 for a given order n based on
% eigenfunctions K and eigenvalues s_\mu
%
% Syntax:  c = outputMatrix_p(n, index, ftm, pickup_pos)
%
% Inputs:
%    n - current order n of Bessel functions
%    index - index \mu for current order n
%    ftm - stuct of tfm parameters 
%    pickup_pos - pickup position
%
% Outputs:
%    c - output matrix for the SSD
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

function c = outputMatrix_p(n, index, ftm, pickup_pos)

c = zeros(1,size(index,1)); 

r = pickup_pos(1);
phi = pickup_pos(2);
z = pickup_pos(3);

nmu = ftm.nmu{n};

for mu = 1:length(c)
    n_real = ftm.n(n); 
    m = index(mu,1);
    nu = index(mu,2); 
    l = ftm.lambda(nu); 
    k = ftm.k(n,m); 
    
    c(mu) = ftm.kprim1(n_real, k, l, r, phi, z).*nmu(mu);
end


end