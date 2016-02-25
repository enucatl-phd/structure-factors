function S = StickyHardSphereStructureFactor(E,p2,L,R,fp,tau,delta)
format long
%%Function to calculate sticky hard sphere structure factor
%as a function of (in order):
%radius, energy, stickiness parameter and volume fraction
%https://kur.web.psi.ch/sans1/sasfit/sasfit.pdf, page 200

% %Parameters to test function
% E = [37.5 34.0 31.0 28.0 25.0 22.0];
% p2 = 2.0E-6;
% L = 0.12;
% R = (0.2E-6:0.2E-6:10.0E-6);
% fp = (0.05:0.05:0.4); %volume fraction
% tau = (0.1:10); %stickiness parameter
% delta = 0; %potential well width

%% Grating interferometer system parameters
lambda = 12.4./E*1e-10; %wavelength [m] for E = 25 keV
theta = p2/L; %theta value for beamline GI

%% Calculating SAXS parameters (Q)
Q = (2*pi./lambda*(theta));
Q = repmat(Q, [length(R) 1 length(tau) length(fp)]);

%% Repmat and permute matrices
tau = repmat(tau, [length(Q(:,1,1,1)) 1 length(Q(1,:,1,1)) length(Q(1,1,1,:))]);
tau = permute(tau, [1 3 2 4]);

R = R.';
R = repmat(R, [1 length(Q(1,:,1,1)) length(Q(1,1,:,1)) length(Q(1,1,1,:))]);
d = R.*2; %diameter

fp = repmat(fp, [length(Q(:,1,1,1)) 1 length(Q(1,1,:,1)) length(Q(1,:,1,1))]);
fp = permute(fp, [1 4 3 2]);

%% Calculating SAXS parameters (S)
A = 2.*Q.*R;

etta = fp.*(((d+delta)./(d)).^3);
epsilon = tau+(etta./(1-etta));
gamma = fp.*((1+(etta./2))./(3.*((1-etta).^2)));
l = (6./etta).*(epsilon-(sqrt((epsilon.^2)-gamma)));
mu = l.*etta.*(1-etta);

alpha = ((1+(2.*etta)-mu).^2)./((1-etta).^4);
beta = -(((3.*etta.*((2+etta).^2))-(2*mu.*(1+(7*etta)+(etta.^2)))...
       +((mu.^2).*(2+etta)))./(2.*((1-etta).^4)));

C = (2.*etta.*l.*sin(A)./A)...
    -(2.*(etta.^2).*(l.^2).*(1-cos(A))./(A.^2))...
    -(((alpha.*(A.^3)).*(sin(A)-(A.*cos(A))))...
    +((beta.*(A.^2)).*((2.*A.*sin(A))-(((A.^2)-2).*cos(A))-2))...
    +(etta.*alpha./2).*((((4.*(A.^3))-(24*A)).*sin(A))...
    -(((A.^4)-(12.*(A.^2))+24).*cos(A))+24))...
    .*(24.*etta./(A.^6));

S = 1./(1-C);