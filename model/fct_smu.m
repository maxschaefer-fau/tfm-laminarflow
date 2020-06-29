%% fct_smu(ftm, index)
% Function to calculate the eigenvalues s_mu
% 
% input:    'ftm'       array of ftm parameters 
%
% output    'smu'       vector of s_mu values 
%
%
%   Maximilian Schäfer @FAU, April 2020

function smu = fct_smu(ftm, diff, index, n)
    smu = zeros(1,size(index,1));
    
    for mu = 1:length(smu)
        m = index(mu,1); 
        nu = index(mu,2);
        
        smu(mu) = -diff.D*(ftm.k(n,m)^2 + ftm.lambda(nu)^2);
    end
end