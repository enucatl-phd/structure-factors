function [DFEC,VR]=DFEC_sphere(D,lambda,corr_length,f_volume,total_thickness,chi)
% function [DFEC,VR]=DFEC(D,lambda)
% Function to calculate the DFEC Coefficient and the Visibility Reduction
% For a spherical model from Appl Opt. 2011 Aug 1; 50(22):
% 4310–4319
%
% D : vector with the diameter of the spheres
% lambda : wavelength
% corr_length : correlation length of the set-up
% f_volume : partial volume of scatters 
% total_thickness : total thickness of the sample
% chi: optical contrast (n_scattering-n_background)
    
n_mic=size(D,2);

for ii=1:n_mic
    Dprime=D(ii)/corr_length;
    if(Dprime>1)
        geo_factor=Dprime_Yashiro(Dprime);
    else
        geo_factor=Dprime;
    end
    prefactor=3*pi^2/lambda^2*f_volume*abs(chi)^2*corr_length;
    DFEC(ii)=prefactor*geo_factor;
    VR(ii)=exp(-DFEC(ii)*total_thickness);
end

end
