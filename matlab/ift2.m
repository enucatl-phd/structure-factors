function g = ift2( G, delta_f )
%ift2 Inverse 2D Fourier transformation of G 
%   Note that equal frequency sampling is assumed for x and y.
N=size(G,1);
g=fftshift(ifft2(ifftshift(G)))*(N*delta_f)^2;
end

