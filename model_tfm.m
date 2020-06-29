%% model_tfm - Derive components of the state-space description
%
% Syntax:  [ftm, state] = model_tfm(diff, ftm)
%
% Inputs:
%    diff - struct containing physical parameters
%    ftm - strcut containing eigenvalues of the tfm
%
% Outputs:
%    ftm - updated struct containing all eigenvalues and eigenfunctions 
%       ftm.smu - eigenvalues s_\mu
%       ftm.nmu - scaling factors 1/N_\mu
%       ftm.k(..) - anonymous functions for eigenfunctions K and \tilde{K}
%    state - struct containing all matrices of the state-space description
%       state.As - open loop state matrix A in the s-domain
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

function [ftm, state] = model_tfm(diff, ftm)
%% Eigenvalues 
ftm.smu = {};

for n = 1:length(ftm.n) 
    ftm.smu{end+1} = fct_smu(ftm, diff, ftm.indexC{n}, n);
end
%% State matrix in s-domain
state.As = {};

for n = 1:length(ftm.n) 
   state.As{end+1} = diag(ftm.smu{n}); 
end

%% scaling factor nmu, i.e. 1/nmu is calculated
ftm.nmu = {};

for n = 1:length(ftm.n)
   ftm.nmu{end+1} = fct_nmu(ftm, diff, ftm.indexC{n},n);
end

%% Anonymous functions for eigenfunctions 
[ftm.kprim1, ftm.kprim1_pm, ftm.kprim2, ftm.kprim3, ftm.kprim4, ftm.kadj1, ftm.kadj2, ...
    ftm.kadj3, ftm.kadj4] = eigenfunctions_cylinder();

end