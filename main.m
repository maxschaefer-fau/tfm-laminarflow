%% main - Main function 
% main function for the numerical evaluation of the proposed model for
% advection-diffusion in a cylinder with laminar flow based on: 
% 
% M. Schäfer, W. Wicke, L. Brand, R. Rabenstein and R. Schober, ?Transfer 
%   Function Models for Diffusion and Laminar Flow in Cylindrical Systems?, 
%   submitted to IEEE Trans. Molecular, Biological and Multi-scale Commun., 
%   2020, [online]: ?
%
% Further information regarding the implementation can be found on: 
%   maximilianschaefer.org/publication/tfm-laminar
%
%
% The code will be updated regularly 
%
% Syntax:  main
%
% Inputs:
%  
% Outputs:
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
% May 2020; Last revision: 28. June 2020

%% Basic definitions
R0 = 1e-4;                                         % Radius of cylinder in [m]
Z0 = 10e-4;                                        % Length of cylinder in [m]
v0 = 0.5e-4;                                       % Flow velocity in [m/s]
D = 12.5e-13;                                          % Diffusion coefficient in [m^2/s]

[diff] = cylinder_parameters(R0, Z0, D, v0);

N = 1;                                            % Number of Bessel functions N, n = 0, ..., N 
M = 30;                                            % Number of roots M, M = 0, ..., M-1
L = 200;                                           % Number of wavenumbers L, \nu = 0, ..., L-1

tdur = diff.Z0_/diff.v0_;                          % simulation duration 
tdur_norm = (tdur/diff.tau);                       % normalized simulation duration
Fs = 5e3;                                          % sampling frequency 

[ftm] = model_parameters(diff, N, M, L);
[ftm, state] = model_tfm(diff, ftm);
[sim] = simulation_parameters(Fs, tdur_norm);

%% Switch between 'diffusion', 'plugflow', 'laminarflow'
% 'diffusion': Pure diffusion with D in the cylinder
% 'plugflow': Diffusion with D and plug flow with v0 in the cylinder
% 'laminarflow': Diffusion with D and laminar flow with v0(1-r^2/R_0^2) in
% the cylinder
 
simulation_state = 'laminarflow';

switch simulation_state
    case 'diffusion' 
        state.As_vr = {};
        for n = 1:length(ftm.n)
           state.As_vr{end+1} = state.As{n}; 
        end
    case 'plugflow'
        state.As_vr = {};
        state.uni = fb_uni(1, ftm.indexC{1}, ftm, diff);
        for n = 1:length(ftm.n)
            state.As_vr{n} = state.As{n} + diff.v0.*state.uni;
        end
    case 'laminarflow'
        state.As_vr = {};
        state.lam = {};
        
        [state.uni, uPhi, uZ] = fb_uni(1, ftm.indexC{1}, ftm, diff);
        
        for n = 1:length(ftm.n)
            state.lam{n} = fb_lam(n, ftm.indexC{n}, ftm, diff, uPhi, uZ, 0);
            state.As_vr{n} = state.As{n} + diff.v0.*state.uni - diff.v0*state.lam{n};
        end
    otherwise
        warning('No simulation state was set')
        return;
end

% clear unused cell arrays to save memory 
clear state.uni state.lam

%% Block decomposition of state matrix 

L = cell(1,ftm.N);
U = cell(1,ftm.N);
Uinv = cell(1,ftm.N);
for n = 1:ftm.N
  [L{n}, U{n}, Uinv{n}] = blockDecomp(state.As_vr{n}, sim.T, 0);  
end
 
% clear unused cell arrays to save memory
clear state.As_vr 

%% Initial conditions yi 
% -- Choose between different scenarios --
% 'zero'
% 'impulseSingle': A Single impulse at a certain position
% 'pointSingle': A Single point at a certain position 
% 'circularSingle': A single circular disk - uniform for r0 = R_0

% initial_state = 'pointSingle'; 
% excite.pos = [0.5, pi/2, 5, 0.4, pi/4, 0.4];       % [r_e, phi_e, z_e, r0, phi0, z0]

initial_state = 'circularSingle';
excite.pos = [0, pi/2, 1, 1*diff.R0, 0.3];        % [r_e, phi_e, z_e, r0, z0]

% initial_state = 'impulseSingle';
% excite.pos = [0, pi/2, 4];        % [r_e, phi_e, z_e]

[excite.yi, excite.norm_yi] = spatialPulse(initial_state, ftm, excite.pos);


%% ---- UNDER CONSTRUCTION ---- %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitation Function fe_x - Spatial pulse shaping

exc_state = 'zero';
excite.pos_fe = [0.6, pi/6, 1, 0.4, pi/4, 0.3];
[excite.fe_x, excite.norm_fe] = spatialPulse(exc_state, ftm, excite.pos_fe);

% Exciation function fe_t - Temporal pulse shaping 

% 'zero': Default setting 
% 'dirac': Instantaneous release at time t_0 
% 'rect': Rectangular release starting at time t_0, width of t_r
% 'raised cosine': Raised cosine centered at time t0, width of tr
excite_state_temp = 'zero';
excite.fe_t = temporalPulse(excite_state_temp, diff, ftm, sim);

% Excitation fe = fe_x*fe_t
excite.fe = cell(1,length(ftm.n));
for n = 1:length(ftm.n)
    excite.fe{n} = excite.fe_t.*excite.fe_x{n};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Numerical evaluation of state equation 
% Exploits the block partioning and uses the block decomposed structure
ybar = cell(1,ftm.N);

for n = 1:ftm.N
[ybar{n}] = simulateTimeDomain(L{n},excite.fe{n},excite.yi{n},...
    sim.T, U{n}, Uinv{n});
end

% clear unused cell arrays to save memory
clear excite.yi excite.fe excite.fe_x excite.fe_t

%% Numerical evaluation of output Equation

pickup.pos = [0, pi/2, 2];
pout = output(ybar, ftm, pickup, sim);

% give current regime 
ad = dimensional_analysis(diff, pickup.pos(3), excite.pos(3));

% normalization factor. 
% Also dependening on receiver from particle based simulations 
% Here for cuboid V_cube = 0.04^3
pickup.norm = 0.04^3;
normalize = (excite.norm_yi{1})/pickup.norm*sim.T; 

pout = real(pout)/normalize;

figure
plot(sim.t*diff.tau, pout); grid on; hold on; 
return 
%% Spatial - r,z 
% Derive a 3D output matrix for 2D plotting over y-z plane
R = 50; 
PHI = pi/2;
Z = 150;
C_rz = cell(1,ftm.N);
for n = 1:ftm.N
    [C_rz{n}, r_rz, phi_rz, z_rz] = spatialMatrix_rz(n, ftm.indexC{n}, ftm, ...
        diff, R, PHI, Z);
end
figure; plot2D_rz(r_rz, z_rz, C_rz, ybar,1,ftm)
animate_rz(r_rz, z_rz, C_rz, ybar, 10, ftm, sim)
