%Dimensional Analysis
%   Source:
%       - V. Jamali, A. Ahmadzadeh, W. Wicke, A. Noel and R. Schober, "Channel Modeling for Diffusive Molecular Communication—A Tutorial Review," in Proceedings of the IEEE, vol. 107, no. 7, pp. 1256-1301, July 2019.
%

function ad = dimensional_analysis(diff, zr, zc)

%%Reference variables
lref = diff.lambda;    %reference length in meter
tref = diff.tau;       %reference duration second

%%System parameters (dimensionless)
D = diff.D;
a = diff.R0;
v0 = diff.v0;
%zc = 0;         %TX center
%zr = 5;         %RX center
d = zr - zc;    %TX-RX distance

%%System parameters (dimensional)
D_ = D*lref^2/tref;
a_ = a*lref;
v0_ = v0*lref/tref;
d_ = d*lref;

%%Regime
dc = a;
veff = v0/2;    %cross-section area-averaged effective flow velocity
%ad = D*d/(veff*a^2);    %Dispersion parameter alpha_d, see tutorial paper
ad = D*d/(veff*dc^2);    %Dispersion parameter alpha_d, see tutorial paper

fprintf('===========================================================================\n')
fprintf('Parameters (dimensionless): D=%.2g, a=%.2g, v0=%.2g, d=%.2g\n', D, a, v0, d)
fprintf('Parameters (dimensional): D=%.2gm^2/s, a=%.2gm, v0=%.2gm/s, d=%.2gm\n', D_, a_, v0_, d_)
fprintf('Dispersion parameter (linear): %.2g (small: flow, large: diffusion)\n', ad)
fprintf('Dispersion parameter (logarithmic): %.2g (negative: flow, positive: diffusion)\n', log10(ad))
fprintf('===========================================================================\n')
end
