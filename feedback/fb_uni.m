%% fb_uni - Derivation of feedback matrix for uniform flow K_uni
%
% Syntax:  [uni, uPhi, uZ] = fb_uni(n, index, ftm, diff)
%
% Inputs:
%    n - current order of Bessel functions n
%    index - index tupels \mu for the order n 
%    ftm - struct containing tfm parameters
%    diff - strcut containing physical parameters
%
% Outputs:
%    uni - feedback matrix K_uni
%    uPhi - part of K_uni containing compenents in phi-direction
%    uZ - part of K_uni containing compenents in z-direction
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

function [uni, uPhi, uZ] = fb_uni(n, index, ftm, diff)

    uni = zeros(ftm.M*ftm.Nu,ftm.M*ftm.Nu);
    nmu = ftm.nmu{n};
    
    % Phi comonent
    uPhi = 2*pi*ones(ftm.M*ftm.Nu);
    
    % z-component
    uZ = zeros(ftm.Nu);
    for nu = 1:length(ftm.nu) 
       for nu_til = 1:length(ftm.nu) 
           
           if(ftm.nu(nu) == 0 && ftm.nu(nu_til) == 0)
                f_z = 0;  
            elseif(ftm.nu(nu_til) == 0 && ftm.nu(nu) ~= 0)
                f_z = 0;
            elseif(ftm.nu(nu) == 0 && ftm.nu(nu_til) ~= 0)
                f_z = -(1/ftm.lambda(nu_til))*...
                    ((-1)^(ftm.nu(nu_til)) - 1);
            elseif(ftm.nu(nu) == ftm.nu(nu_til))
                f_z = 0;
            else
                f_z = - ftm.lambda(nu_til)/(ftm.lambda(nu_til)^2 - ftm.lambda(nu)^2)*...
                    ((-1)^(ftm.nu(nu_til) + ftm.nu(nu)) - 1);
           end
           
           uZ(nu,nu_til) = f_z; 
       end
    end
    
    % r-component 
    uR = zeros(ftm.M*ftm.Nu); 
    
    for m = 1:length(ftm.m) 
        for m_til = 1:length(ftm.m) 
            if(m ~= m_til) 
                f_r = 0;
            elseif(ftm.k(n,m) == 0 && ftm.n(n) >= 1)
                f_r = 0;
            elseif(ftm.k(n,m) == 0 && ftm.n(n) == 0)
                f_r = diff.R0^2/2;
            else
                f_r = diff.R0^2/2*(1 - ftm.n(n)^2/ftm.k(n,m)^2)*besselj(ftm.n(n),ftm.k(n,m))^2; 
            end
            uR(m,m_til) = f_r; 
        end 
    end
    
    for mu = 1:size(uni,1) 
        m = index(mu,1); 
        nu = index(mu,2);      
       for mu_til = 1:size(uni,1) 
            m_til = index(mu_til,1); 
            nu_til = index(mu_til,2);
            
            uni(mu,mu_til) = uPhi(mu,mu_til)*uZ(nu,nu_til)*...
                uR(m,m_til)*nmu(mu_til)*ftm.lambda(nu); 
       end
    end
end