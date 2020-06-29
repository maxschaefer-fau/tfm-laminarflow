function nmu = fct_nmu(ftm, diff, index, n)

nmu = zeros(1,size(index,1)); 

for mu = 1:length(nmu)
   m = index(mu,1); 
   nu = index(mu,2); 
   
   % phi component
   nmu_phi = 2*pi; 
   nmu_phi_inv = 1/nmu_phi; 
   
%    z-component 
   if(ftm.lambda(nu) == 0)
      nmu_z = 0; 
      nmu_z_inv = 0; 
   else
      nmu_z = diff.Z0/2; 
      nmu_z_inv = inv(nmu_z); 
   end

   % r-component
   if(ftm.k(n,m) == 0)
       nmu_r_inv = 2/diff.R0;
   else
       nmu_r_inv = 2/diff.R0*ftm.k(n,m)^2/(...
        besselj(ftm.n(n),ftm.k(n,m))^2* ... 
        (ftm.k(n,m)^2 - ftm.n(n)^2/diff.R0^2)...
        );
   end
   
   nmu(mu) = nmu_phi_inv*nmu_z_inv*nmu_r_inv;
end


end