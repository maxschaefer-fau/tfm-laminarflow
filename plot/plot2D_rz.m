%% plot2D_rz - plots the 2D r-z concentration at time t0
%
% 
% Syntax:  plot2D_rz(r, z, space, time, t0, ftm)
%
% Inputs:
%    r - sampling points in r-direction  
%    z - sampling points in z-direction  
%    space - array of eigenfunctions [x,y,modes]
%    time - system states [time, modes]
%    t0 - time position for the plot
%    ftm - struct with tfm parameters
%
% Outputs:
%       
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
% Authors: Dr.-Ing. Maximilian Schaefer, Dr.-Ing. Sebastian Schlecht 
% University of Erlangen-Nuremberg
% email address: max.schaefer@fau.de
% Website: maximilianschaefer.org
% 28. June 2020; Last revision: 28. June 2020

function plot2D_rz(r, z, space, time, t0, ftm)

time = permute(time,[3 2 1]);

h = surf(z, r, zeros(length(r),length(z)),'edgecolor','none');
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'fontsize',18)
zlim([0 1e9]);
xlabel('Space $[z]$','interpreter','latex');
ylabel('Space $[y]$','interpreter','latex');
yticks([-1 -0.5 0 0.5 1])
yticklabels({'$-R_0$','$-\frac{R_0}{2}$','$0$','$\frac{R_0}{2}$','$R_0$'})
view([0 90]);
set(gca,'DataAspectRatio',[2 1 1])
shading interp;

dat = zeros(length(r),length(z));
for n = 1:length(space)
    time_n = time{n}(:,t0);
    time_n = permute(time_n,[3 2 1]);
    d = sum(space{n}.*time_n, 3);
    % Consider negative terms for n \neq 0
    if(ftm.n(n) ~= 0)
        dat = dat + 2*real(d);
    else
        dat = dat + d;
    end
end
set(h, 'ZData', real(dat));


