function S = HardSphereStructureFactor2(Q,R,f)

%%Function to calculate structure factors for different radii hard spheres
%D.J. Kinning et al., Macromolecules 17 (1984) 1712

% %Parameters to test function
% Q = momentum;
% R = radius of spheres;
% f = fraction volume;


%% Calculating SAXS parameters (Q, S)

A = 2*Q*R;

alpha = ((1+(2.*f)).^2)./((1-f).^4);
beta = -6.*f.*((1+(f./2)).^2)./((1-f).^4);
gamma = ((0.5.*(f)).*((1+(2.*f)).^2))./((1-f).^4);
            
G1 = (alpha./(A.^2)) .* (sin(A)-(A.*cos(A)));
G2 = (beta./(A.^3)) .* ((2.*A.*sin(A))+((2-(A.^2)).*cos(A))-2);
G3 = (gamma./(A.^5)) .* (((-A.^4).*cos(A)) + (4.*(((3*(A.^2)-6).*cos(A))+((sin(A).*((A.^3)-(6.*A))))+6)));
G = (G1+G2+G3);
            
S = 1./(1+(24.*f.*(G./A))); %as a function of: (in order)
                            %radius,energy,volume fraction

[n,m]=size(A);
S(A==0)=1/(1+24*f*(alpha/3+beta/4+gamma/6));

end








