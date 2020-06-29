%% fct_smuIndex(ftm)
% Function to create index-tupels 'mu' = ['n','m','nu'] as 
% s_mu = -(k_nm)^2 - (lambda_nu)^2
%
% From now on the index 'mu' can be used instead of all three indices. 
% (e.g. one for-loop instead of three, more clear code)
%
% input:    'ftm'       array of ftm parameters 
%
% output    'index'    vector of mu = 1:Mu index tupels 
%
%
%   Maximilian Schäfer @FAU, April 2020

function index = fct_index(ftm)

index = zeros(ftm.Mu,3);

% nu index
nu = 1:ftm.Nu;
index(:,3) = repmat(nu,1,ftm.N*ftm.M);

% m index 
index(:,2) = repmat(repelem(ftm.m, ftm.Nu),1,ftm.N);

% n index
n = 1:ftm.N;
index(:,1) = repelem(n,ftm.M*ftm.Nu);
end