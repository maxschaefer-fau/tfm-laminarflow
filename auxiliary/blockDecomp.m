%% blockDecomp - Block decomposition of closed loop state matrix
% Closed loop state matrix is decomposed and transformed into the discrete
% time domain
%
% Syntax:  [L, U, Uinv] = blockDecomp(As, T, saveMat)
%
% Inputs:
%    As - closed loop state matrix in the s-domain
%    T - sampling time
%    saveMat - flag wether the matrices should be saved (not recommended)
%
% Outputs:
%    L - diagonal matrix of eigenvalues in the discrete-time domain
%    U - matrix of eigenvectors 
%    Uinv - inverse matrix of eigenvectors 
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

function [L, U, Uinv] = blockDecomp(As, T, saveMat)

    Lblock = zeros(size(As));
    Ublock = zeros(size(As));
    Uinvblock = zeros(size(As));
    
    [U,Ls] = eig(As);
    L = expm(Ls*T);
    Uinv = inv(U);

    if(saveMat == 1)
       save('Decomp.mat','Ublock', 'Uinvblock','Lblock','-v7.3');
    end
end