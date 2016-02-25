function sigma =Sigma_VR_detector(Nps,a0r,vr,f1_r,f1_s,T,Vratio)
% function sigma =Sigma_VR(Nps,a0r,vr,f1_r,f1_s,T,Vratio)
% Function to calculate the visibility reduction error
% from Rev. Sci. Instrum. 81, 073709 (2010)
%
% Nps   : number of phase steps
% a0r   : intensity for a single phase step without sample
% vr    : Visibility reduction without sample (a1r)
% f1_r  : signal to noise transfer Eq. 7 without sample
% f1_s  : signal to noise transfer Eq. 7 with sample
% T     : Transmission ratio after the sample
% Vratio: Visibility ratio due to the sample

I_tot=Nps*a0r;
sigma=f1_r/vr^2/I_tot*(vr^2*(1+f1_s/f1_r./T)+2*(1+f1_s/ ...
                                                  f1_r./T./ ...
                                                Vratio.^2));
sigma=sqrt(Vratio.^2.*sigma);
end
