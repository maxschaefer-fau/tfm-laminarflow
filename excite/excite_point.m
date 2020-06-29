%% excite_point - initial values \bar{y}_init for a point distribution
%
% Syntax:  [excite, norm] = excite_point(ftm, index, n, excite_pos)
%
% Inputs:
%    ftm - struct containing all eigenvalues and eigenfunctions of tfm
%    index - index \mu for given order n 
%    n - current order n of Bessel functions
%    excite_pos - excitation position in the cylinder
%
% Outputs:
%    excite - initial values for the point distribution
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

function [excite, norm] = excite_point(ftm, index, n, excite_pos)

re = excite_pos(1); 
phie = excite_pos(2); 
ze = excite_pos(3); 
r0 = excite_pos(4);
phi0 = excite_pos(5);
z0 = excite_pos(6);
 
kz = 2*pi/z0;
kr = 2*pi/r0;
kphi = 2*pi/phi0;

% Functions for integration
fun_r = @(r,n,k,kr,re) besselj(n,k*r).*r.*0.5.*(1 + ...
        cos(kr*(r - re)));
fun_z = @(z,l,kz,ze)sin(l*z).*0.5.*(1 + ...
        cos(kz*(z - ze)));
    
fun_phi = @(phi, n, phie, kphi) exp(-1j*n*phi).*0.5.*(1 + ...
        cos(kphi*(phi - phie))); 
    

excite = zeros(size(index,1),1); 

f_phi = integral(@(phi)fun_phi(phi, ftm.n(n), phie, kphi), phie-phi0/2, phie+phi0/2);

f_r = zeros(1,length(ftm.m));
for m = 1:length(ftm.m)
    k = ftm.k(n,m);
    f_r(m) = integral(@(r)fun_r(r, ftm.n(n), k, kr,re), re-r0/2, re+r0/2);
end

f_z = zeros(1,length(ftm.nu));
for nu = 1:length(ftm.nu)
    l = ftm.lambda(nu);
    f_z(nu) = integral(@(z)fun_z(z,l,kz,ze), ze-z0/2, ze+z0/2);
end


for mu = 1:length(excite)
    m = index(mu,1);
    nu = index(mu,2); 

    excite(mu) = f_r(m)*f_phi*f_z(nu);
end

% Normalization 
funr = @(r,re,kr) r.*0.5.*(1 + cos(kr*(r - re)));
funz = @(z,kz,ze) 0.5.*(1 + cos(kz*(z - ze)));
funphi = @(phi,kphi,phie) 0.5.*(1 + cos(kphi*(phi - phie)));

nphi = integral(@(phi)funphi(phi,kphi,phie),phie-phi0/2, phie+phi0/2);
nr = integral(@(r)funr(r,re,kr),re-r0/2, re+r0/2);
nz = integral(@(z)funz(z,kz,ze),ze-z0/2, ze+z0/2);

norm = nr*nz*nphi;
end

