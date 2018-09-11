% fresnel.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs the Fresnel transform reconstruction
% and outputs the electric field strength as a matrix of complex numbers b.

function [b] = fresnel(hol,d,lambda,res,ref)

% Declaring pixel and sensor size in the hologram plane:
L = size(hol);
M = L(2);      % number of pixels along x-direction
N = L(1);      % number of pixels in y-direction
del_xi = res;  % pixel size in x-direction
del_eta = res; % pixel size in y-direction
chirp = ones(N,M);

% Computing the length scales in the reconstruction plane:
xsize_fres = abs(d*lambda/M/res);
ysize_fres = abs(d*lambda/N/res);
xfres = [-(M/2)*xsize_fres:xsize_fres:(M/2)*xsize_fres];
xfres = xfres(1:M);
yfres = ([-(N/2)*ysize_fres:ysize_fres:(N/2)*ysize_fres]);
yfres = yfres(1:N);

k = 2*pi/lambda;

% Generating the chirp function required by the Fresnel transform:
status = sprintf('Generating chirp...')
for n = 0:N-1
        for m = 0:M-1
            
            % Chirp function:
            chirp(n+1,m+1) = exp((1i*pi/(d*lambda))*(((n-N/2)^2*...
                del_eta^2)+((m-M/2)^2*del_xi^2)));

            % Complex, constant phase factor:
            % This factor is only required if we want a correctly scaled
            % reconstructed intensity distribution.  It is not needed in
            % computing the interference phase distribution.
            % A(n+1,m+1) = (exp(1i*k*d)/(1i*lambda*d))*...
            % exp((1i*pi/(d*lambda))*(xfres(m+1)^2+yfres(n+1)^2));
            
        end
end

status = sprintf('done with chirp. \n')

% Set the complex factor to 1 when only reconstructing phase.
A = 1; 

% Multiply the hologram intensity distribution with the reference beam:
hol_ref = hol*ref;

% Convolving the product of the hologram function and reference wave with 
% the chirp function.
status = sprintf('Computing Fresnel transform...')
b = A.*ifftshift(ifft2(fftshift(hol_ref.*chirp)));
status = sprintf('done with Fresnel transform. \n')