%% cylinder_parameters - Definition of geometrical cylinder parameters
%Definition of radius R, and length Z0 of the cylinder
%Definition of diffusion coefficient D, and flow velocity v0 = v_eff/s
%
% Syntax:  [diff] = cylinder_parameters(R0, Z0, D, v0)
%
% Inputs:
%    R0 - Radius of the cylinder
%    Z0 - Length of the cylinder considered for numerical evaluation
%    D - Diffusion coefficient
%    v0 - Flow velocity 
%
% Outputs:
%    diff - struct containing all physical parameters of the cylinder
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

function [diff] = cylinder_parameters(R0, Z0, D, v0)

% Geometrical parameters
diff.R0_ = R0;                                         % Radius of cylinder in [m]
diff.Z0_ = Z0;                                        % Length of cylinder in [m]

% Diffusion and flow
diff.v0_ = v0;                                       % Flow velocity in [m/s]
diff.D_ = D;                                          % Diffusion coefficient in [m^2/s]

% Normalization - Reference length, Reference time 
diff.tau = 1e2;                                          % Reference time in [s]
diff.lambda = diff.R0_;                                  % Reference length in [m]

% Apply normalization
diff.R0 = diff.R0_/diff.lambda;                          % Normalized cylinder radius
diff.Z0 = diff.Z0_/diff.lambda;                          % Normalized cylinder length
diff.v0 = diff.v0_/diff.lambda*diff.tau;                 % Normalized flow velocity
diff.D = diff.D_/diff.lambda^2*diff.tau;                 % Normalized diffusion coefficient


end