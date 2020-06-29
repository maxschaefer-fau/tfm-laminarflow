%% excite_circular - initial values \bar{y}_init for a cicular distribution 
%  is a uniform distribution for r0 = R0
%
% Syntax:  [excite, norm] = excite_circular(ftm, index, n, excite_pos)
%
% Inputs:
%    ftm - struct containing all eigenvalues and eigenfunctions of tfm
%    index - index \mu for given order n 
%    n - current order n of Bessel functions
%    excite_pos - excitation position in the cylinder
%
% Outputs:
%    excite - initial values for the circular distribution
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

function [excite, norm] = excite_circular(ftm, index, n, excite_pos)

ze = excite_pos(3); 
r0 = excite_pos(4);
z0 = excite_pos(5);

k0 = 2*pi/excite_pos(5); 


excite = zeros(size(index,1),1); 

% Functions for integration
fun_r = @(r,n,k) besselj(n,k*r).*r;
fun_phi = @(phi,n) exp(-1j*n*phi);
fun_z = @(z,l,k0,ze)sin(l*z).*0.5.*(1 + ...
        cos(k0*(z - ze)));

f_phi = integral(@(phi)fun_phi(phi, ftm.n(n)), -pi, pi);
    
f_r = zeros(1,length(ftm.m));
for m = 1:length(ftm.m)
    k = ftm.k(n,m);
    f_r(m) = integral(@(r)fun_r(r, ftm.n(n), k), 0, r0);
end

f_z = zeros(1,length(ftm.nu));
for nu = 1:length(ftm.nu)
    l = ftm.lambda(nu);
    f_z(nu) = integral(@(z)fun_z(z,l,k0,ze), ze-z0/2, ze+z0/2);
end

for mu = 1:length(excite) 
    
    m = index(mu,1);
    nu = index(mu,2); 
   
    excite(mu) = f_r(m)*f_phi*f_z(nu);
end

% normalization
funr = @(r) r;
funz = @(z,ze)0.5.*(1 + cos(k0*(z - ze)));

nr = integral(@(r)funr(r), 0, r0);
nphi = 2*pi;
nz = integral(@(z)funz(z, ze),ze-z0/2, ze+z0/2);

norm = nr*nz*nphi;
    
end