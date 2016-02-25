function Pfactor = SphereFormFactor2(Q,R)
%%Function to calculate the sphere form factor for different radii hard spheres
%D.J. Kinning et al., Macromolecules 17 (1984) 1712

% Q=4*pi/lambda*sin(theta/2)
% R= radius of microspheres

A = Q*R;
Pfactor = (3./(A.^3)).*(sin(A)-(A.*cos(A))); 
Pfactor(A==0)=1/3;
end
