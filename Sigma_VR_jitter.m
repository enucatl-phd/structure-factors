function sigma =Sigma_VR_jitter(Nps,p2,sigma_x,vr,Vratio)
% function sigma =Sigma_VR_jitter(Nps,p2,vr,Vratio)
% Function to calculate the visibility reduction error due to jitter
% from Rev. Sci. Instrum. 81, 073709 (2010)
%
% Nps     : number of phase steps
% p2      : period of g2
% sigma_x : error positioning
% vr      : Visibility reduction without sample (a1r)
% Vratio  : Visibility ratio due to the sample

sigma=2*pi^2/Nps*(sigma_x/p2)^2*(2+vr^2*(1+Vratio));
sigma=sqrt(Vratio.^2.*sigma);
end
