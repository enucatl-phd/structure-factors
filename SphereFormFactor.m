function F = SphereFormFactor(E,p2,L,R)

%%Function to calculate the sphere form factor for different radii hard spheres
%D.J. Kinning et al., Macromolecules 17 (1984) 1712

% %Parameters to test function
% E = [37.5 34.0 31.0 28.0 25.0 22.0];
% p2 = 2.0E-6;
% L = 0.12;
% R = [5.00E-6 3.2E-6 1E-6];


%% Grating interferometer system parameters
lambda=12.4./E*1e-10; %wavelength for E = 25 keV [m]
theta = p2/L; %theta value for beamline GI


%% Sample parameters
%d = R.*2;


%% Calculating SAXS parameters (Q, S)

Q = (2*pi./lambda*(theta));
Q = repmat(Q, [length(R) 1]);
R = repmat(R, [length(Q(1,:)) 1]).';

A = Q.*R;
            
F = (3./(A.^3)).*(sin(A)-(A.*cos(A))); %as a function of in order: radius,energy,volume fraction