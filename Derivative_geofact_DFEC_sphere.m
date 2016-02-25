function factor=Derivative_geofact_DFEC_sphere(D)
% function factor=Geofact_DFEC_sphere(D)
% Function to calculate the derivative of the geometrical factor of the DFEC
% Coefficient according to 
% Appl Opt. 2011 Aug 1; 50(22): 4310–4319
% when D is bigger than 1
% D : Diameter of the spheres/ correlation_length fo the GI

    
A1=sqrt(D^2-1);
A2=(1+1/2/D^2);
A1prime=D/A1;
A2prime=-1/D^3;
B1=(1/D-1/4/D^3);
B2=log((D+A1)/(D-A1));
B1prime=-1/D^2+3/4/D^4;
B2prime=(D-A1)/(D+A1)*((1+A1prime)/(D-A1)-(D+A1)/(D-A1)^2*(1-A1prime));
factor=1-A1prime*A2-A1*A2prime+B1prime*B2+B2prime*B1;
end
