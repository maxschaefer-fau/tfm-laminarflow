%% model_parameters - Definition basic parameters for the transfer function model
% Definition of N, M, L 
% Derivation of wavenumbers \lambda_\nu, and roots k_{n,m}
%
% Syntax:  [ftm] = model_parameters(N, M, L)
%
% Inputs:
%    N - maximum number of Bessel function orders
%    M - number of roots 
%    L - number of wave-numbers \lambda_\nu
%
% Outputs:
%    ftm - struct containing values k_{n,m}, \lambda_\nu
%       ftm.Nu - number of wavenumbers L
%       ftm.M - number of roots M
%       ftm.N - maximum number of Bessel function orders
%       ftm.n - array of Bessel order indices 
%       ftm.m - array of root numbers indices
%       ftm.nu - array of wave-number indices 
%       ftm.lambda - array of wavenumbers \lambda_\nu
%       ftm.k - array of roots k_{n,m}
%       ftm.index - index tupels \mu -> n,m,\nu
%       ftm.indexC - index tupels \mu as cell array - exploits block
%       partioning 
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

function [ftm] = model_parameters(diff, N, M, L)
ftm.Nu = L;                               % number of z-eigenvalues (also named L in paper)
ftm.M = M;                                % number of roots per order
ftm.N = N;                                % maximum order of bessel-functions 

ftm.Mu = ftm.N*ftm.M*ftm.Nu;                % Total number of eigenvalues (also named Q in paper) 

% Arrays of indices 
ftm.nu = 0:ftm.Nu - 1;                      % All indices nu = 0:Nu-1
ftm.m = 1:ftm.M;                            % All indices m = 1:M
ftm.n = 0:ftm.N-1;                            % All indices n = 0:N-1

% Index all eigenvalues -- Not necessary
ftm.mu = 1:ftm.Mu;                          % All indies mu = 1:Mu
%% Eigenvalues s_mu
[ftm.lambda] = fct_lambda(ftm, diff);                                      % Calculate lambda_nu values
[ftm.k, ] = BessDerivZerosBisect2(ftm.n,ftm.m);                            % Calcualte bessel zeros k_nm

% Add roots at zero, i.e., k = 0 for order n = 0
tmp = ftm.k(1,:); 
tmp = [0 tmp(1:end-1)];
ftm.k(1,:) = tmp;

% create the index 'mu' as index-tupel of 'n,m,nu'
ftm.index = fct_index(ftm);                                                % Full index                       
ftm.indexC = {};                                                           % Exploit block partitioning

for n = 1:length(ftm.n)
    foo = zeros(ftm.M*ftm.Nu,2);
    nu = 1:ftm.Nu;
    foo(:,2) = repmat(nu,1,ftm.M);
    foo(:,1) = repelem(ftm.m, ftm.Nu);
    
    ftm.indexC{n} = foo;
end

end