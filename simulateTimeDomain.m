%% simulateTimeDomain - numerical evaluation of modified state equation
% Evaluates the modified state equation based on the block decomposition of
% state matrix Ac - Function is called N times for each individual order n 
%
% Syntax:  [ybar] = simulateTimeDomain(Az,fe,yi,T, U, Uinv)
%
% Inputs:
%    Az - state Ac_n for current order n 
%    fe - external sources 
%    yi - initial values 
%    T - simulation time step
%    U - matrix of eigenvectors from block decomposition
%    Uinv - inverse matrix of eigenvectors from block decomposition
%
% Outputs:
%    ybar - system states ybar_n for current order n
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

function [ybar] = simulateTimeDomain(Az,fe,yi,T, U, Uinv)

Mu = size(Az,1);
uybar = zeros(Mu,size(fe,2));

% extend initial conds:
yi_time = zeros(Mu, size(fe,2));
yi_time(:,1) = yi;

excite = fe + yi_time; 
uexcite = Uinv*excite;

for mu = 1:Mu
    uybar(mu,:) = filter(1,[1, -Az(mu,mu)],T*uexcite(mu,:));
    % O(sim.t)
end
% 1 Addition + 1 Multiplikation pro Zeitschritt und pro Mode
% O(sim.t*Mu)

% u^-1*ybar --> ybar
ybar = U*uybar;
% NxN ist O(N^3) -- (Mu x Mu)*(Mu x sim.t) O(Mu*Mu*sim.t) 