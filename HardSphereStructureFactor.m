function S = StructureFactor(E,p2,L,R,etta)

%%Function to calculate structure factors for different radii hard spheres
%D.J. Kinning et al., Macromolecules 17 (1984) 1712

% %Parameters to test function
% E = [37.5 34.0 31.0 28.0 25.0 22.0];
% p2 = 2.0E-6;
% L = 0.12;
% R = [5.00E-6 3.2E-6 1E-6];
% etta = [0.23 0.25 0.5];


%% Grating interferometer system parameters
lambda=12.4./E*1e-10; %wavelength for E = 25 keV [m]
theta = p2/L; %theta value for beamline GI


%% Sample parameters
%d = R.*2;


%% Calculating SAXS parameters (Q, S)

Q = (2*pi./lambda*(theta));
Q = repmat(Q, [length(R) 1]);
R = repmat(R, [length(Q(1,:)) 1]).';

A = 2.*Q.*R;

alpha = ((1+(2.*etta)).^2)./((1-etta).^4);
beta = -6.*etta.*((1+(etta./2)).^2)./((1-etta).^4);
gamma = ((0.5.*(etta)).*((1+(2.*etta)).^2))./((1-etta).^4);

A = repmat(A, [1 1 length(alpha)]);
Q = repmat(Q, [1 1 length(alpha)]);
R = repmat(R, [1 1 length(alpha)]);


alpha = repmat(alpha, [length(A(:,1,1)) 1 length(A(1,:,1))]);
alpha = permute(alpha, [1 3 2]);
beta = repmat(beta, [length(A(:,1,1)) 1 length(A(1,:,1))]);
beta = permute(beta, [1 3 2]);
gamma = repmat(gamma, [length(A(:,1,1)) 1 length(A(1,:,1))]);
gamma = permute(gamma, [1 3 2]);
etta = repmat(etta, [length(A(:,1,1)) 1 length(A(1,:,1))]);
etta = permute(etta, [1 3 2]);
          
            
G1 = (alpha./(A.^2)) .* (sin(A)-(A.*cos(A)));
G2 = (beta./(A.^3)) .* ((2.*A.*sin(A))+((2-(A.^2)).*cos(A))-2);
G3 = (gamma./(A.^5)) .* (((-A.^4).*cos(A)) + (4.*(((3*(A.^2)-6).*cos(A))+((sin(A).*((A.^3)-(6.*A))))+6)));
G = (G1+G2+G3);
            
S = 1./(1+(24.*etta.*(G./A))); %as a function of: (in order) radius,energy,volume fraction