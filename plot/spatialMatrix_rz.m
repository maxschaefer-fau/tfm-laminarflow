%% spatialMatrix_rz - output matrix for complete r-z domain
% output matrix for the r-z domain to create 2D plots of particle
% distribution
%
% Syntax:  [Cper, r, phi, z] = spatialMatrix_rz(n, index, ftm, diff, R, PHI, Z)
%
% Inputs:
%    n - current order n of Bessel function
%    index - index \mu for current order n 
%    ftm - struct of tfm parameters 
%    diff - struct of physical parameters 
%    R - number of sampling points in r-direction 
%    PHI - phi-position (has to be fixed)
%    Z - number of sampling points in z-direction 
%
% Outputs:
%    Cper - output matrix for the r-z-plot (very large)
%    r - sampling points in r-direction 
%    phi - phi-position
%    z - sampling points in z-direction 
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

function [Cper, r, phi, z] = spatialMatrix_rz(n, index, ftm, diff, R, PHI, Z)

deltaR = diff.R0/R;
r = -diff.R0:deltaR:diff.R0;

phi = PHI;

deltaZ = diff.Z0/Z;
z = 0:deltaZ:diff.Z0;

ri = 1:length(r);
zi = 1:length(z);

Cper = zeros(size(index,1),length(ri), length(zi)); 
k = ftm.k(n,index(:,1)).';
l = ftm.lambda(index(:,2)).';

Cper(:,:,:) = besselj(ftm.n(n),k.*r(ri))...
                .*exp(1j*ftm.n(n).*phi)...
                .*sin(l.*permute(z(zi),[1 3 2]));

Cper = transpose(ftm.nmu{n}).*Cper;
Cper = permute(Cper,[2 3 1]);
end