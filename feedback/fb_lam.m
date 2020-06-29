%% fb_lam - Derivation of feedback matrix for parabolic profile K_lam
%
% Syntax:  [lam] = fb_lam(n, index, ftm, diff, uPhi, uZ, saveMat)
%
% Inputs:
%    n - current order of Bessel functions n
%    index - index tupels \mu for the order n 
%    ftm - struct containing tfm parameters
%    diff - strcut containing physical parameters
%    uPhi - part of K_uni containing compenents in phi-direction (similiar
%    for K_lam)
%    uZ - part of K_uni containing compenents in z-direction (similiar
%    for K_lam)
%    saveMat - flag wether the matrix should be saved (not recommended)
%
% Outputs:
%    lam - feedback matrix K_lam
%       
%
% Other m-files required: none
% Subfunctions: 
%       [f_r_lam] = orderZero(diff, kqi, knm, m, m_til)
%       Function to derive the r-component of K_lam for n = 0 (possible in
%       closed form, no integration necessary)
%
% MAT-files required: none
%
% See also: none
% Author: Dr.-Ing. Maximilian Schaefer, 
% University of Erlangen-Nuremberg
% email address: max.schaefer@fau.de
% Website: maximilianschaefer.org
% 28. June 2020; Last revision: 28. June 2020

function [lam] = fb_lam(n, index, ftm, diff, uPhi, uZ, saveMat)

    lam = zeros(ftm.M*ftm.Nu,ftm.M*ftm.Nu);
    nmu = ftm.nmu{n};
    tr_r = 0:1e-3:diff.R0;

    % Phi and z component already done in fb_uni 
    % saved in uPhi, uZ
    
    uR = zeros(ftm.M, ftm.M); 
    
    for m = 1:length(ftm.m)
       for m_til = 1:length(ftm.m)
           
           % order zero
           if(ftm.n(n) == 0)
              knm = ftm.k(n,m); 
              kqi = ftm.k(n,m_til); 
              uR(m,m_til) = orderZero(diff, kqi, knm, m, m_til);
           else 
              knm = ftm.k(n,m); 
              kqi = ftm.k(n, m_til);
              n_real = ftm.n(n);
              
              tr_fun = besselj(n_real,knm.*tr_r).*besselj(n_real,kqi.*tr_r).*tr_r...
                     .*(tr_r.^2/diff.R0.^2); 
              uR(m,m_til) = trapz(tr_r,tr_fun);
           end
       end
    end
    
    for mu = 1:size(lam,1) 
        m = index(mu,1); 
        nu = index(mu,2);      
       for mu_til = 1:size(lam,1) 
            m_til = index(mu_til,1); 
            nu_til = index(mu_til,2);
            
            lam(mu,mu_til) = uPhi(mu,mu_til)*uZ(nu,nu_til)*...
                uR(m,m_til)*nmu(mu_til)*ftm.lambda(nu); 
       end
    end
    
if(saveMat == 1)
   save('FB_lam.mat','lam','-v7.3');
end
    
end


function [f_r_lam] = orderZero(diff, kqi, knm, m, m_til)
    if(m == m_til && kqi == 0 && knm == 0)
       f_r_lam = 0.25*diff.R0^4; 
    elseif(m == m_til)
       f_r_lam = diff.R0^4/6*besselj(0,knm*diff.R0)^2 + ...
           diff.R0^3/3*besselj(0,knm*diff.R0)*besselj(1,knm*diff.R0) + ...
           (diff.R0^4/6 - diff.R0^2/3)*besselj(1,knm*diff.R0)^2;
    else
       f_r_lam = 2*diff.R0^2/(knm^2 - kqi^2)^2*...
           ((knm^2 + kqi^2)*besselj(0,knm*diff.R0)*besselj(0,kqi*diff.R0) + ...
           2*knm*kqi*besselj(1,knm*diff.R0)*besselj(1,kqi*diff.R0)) + ...
           (4*diff.R0*(knm^2 + kqi^2)/(knm^2 - kqi^2)^3 - diff.R0^3/(knm^2 - kqi^2))*...
           (kqi*besselj(0,knm*diff.R0)*besselj(1,kqi*diff.R0) ...
           - knm*besselj(1,knm*diff.R0)*besselj(0,kqi*diff.R0));
    end
end