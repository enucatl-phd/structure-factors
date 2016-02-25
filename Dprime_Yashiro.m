function factor=Dprime_Yashiro(D)
% function factor=Dprime_Yashiro(D)
% Function to calculate the geometrical factor of the DFEC
% Coefficient according to 
% Appl Opt. 2011 Aug 1; 50(22): 4310–4319
% when D is bigger than 1
% D : Diameter of the spheres/ correlation_length fo the GI

A=sqrt(D^2-1)*(1+1/D^2/2);
B=1/D-1/4/D^3;
C=log((D+sqrt(D^2-1))/(D-sqrt(D^2-1)));

factor=D-A+B*C;
end
