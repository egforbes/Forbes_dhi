dr = 0.25;
C = zeros(10,10);
for j = 1:1:10
  for i = 1:1:j
     C(i,j) = sqrt(((j+1))^2-(i)^2)-sqrt((j)^2-(i)^2);
  end
end
C=C*2*dr;
C

break;


width = 2;
s = zeros(100,100);
sf = boxcar2(s,width); 
size(sf)
phase_f = sf(2*width+1:end-(2*width+1),:);
size(phase_f)

break;

status = sprintf('Initializing... \n')
  
res = (3.85e-6);
lambda = 532e-9;
d =.50; % [m] Reconstruction distance
n_c = 2.312e27; % [m^-3] Critical plasma density
ref = 1;
hol_def_rgb = imread('../../data/Base_Shot_26.JPG');
hol_base_rgb = imread('../../data/Plasma_Shot_28.JPG');

rgb_red = 0.3;
rgb_green = 0.59;
rgb_blue = 0.11;
    
hol_def = double(hol_def_rgb(:,:,1).*rgb_red + hol_def_rgb(:,:,2).*rgb_green + hol_def_rgb(:,:,3).*rgb_blue)';

hol_base = double(hol_base_rgb(:,:,1).*rgb_red + hol_base_rgb(:,:,2).*rgb_green + hol_base_rgb(:,:,3).*rgb_blue)';

del_xi = res;  % pixel size in x-direction
del_eta = res; % pixel size in y-direction

[phase_diff,x,y,xsize_fres,ysize_fres] = DHI_rec(hol_base,hol_def,d,res);

MS_Real = imag(phase_diff);
fidR = fopen("../../data/phaseDiffMatlab.dat",'w');
fwrite(fidR,MS_Real,'double');
fclose(fidR);

image(phase_diff);

break;

M = size(hol_def,2);
N = size(hol_def,1);

x_plt = -max(x)/2:x(2)-x(1):max(x)/2;
y_plt = -max(y)/2:y(2)-y(1):max(y)/2;

% Extracting the twin image:

% Labeling the matrix indices:
x_ind_vec = 1:1:M;
y_ind_vec = 1:1:N;

xmin=0;
xmax=0;
ymin=0;
ymax=0;
 
% Interpolating to match a spatial coordinate values to the
% nearest integer matrix indices:
x_min_ind = round(interp1(x,x_ind_vec,xmin));
x_max_ind = round(interp1(x,x_ind_vec,xmax));
y_min_ind = round(interp1(y,y_ind_vec,ymin));
y_max_ind = round(interp1(y,y_ind_vec,ymax));

% Defining the x,y spatial coordinate vectors for the extracted twin image:
x_twin = x(x_min_ind:x_max_ind);
x_twin = x_twin(1:10:end);
y_twin = y(y_min_ind:y_max_ind);
y_twin = y_twin(1:10:end);

twin_img = phase_diff(y_min_ind:y_max_ind,x_min_ind:x_max_ind);
twin_img = twin_img(1:10:end,1:10:end);

% Smoothing and unwrapping the twin image:
width = 8;
[unwrapped,phase_f,unsmoothed_unwrapped] = smooth_unwrap(twin_img,width);





break;



if (0)

status = sprintf('Computing Fresnel transform of def...')
b_def =fft2(hol_def.*chirp);
b_def =fftshift(b_def);

b_base =fft2(hol_base.*chirp);
b_base =fftshift(b_base);

phase_diff = angle(b_base.*conj(b_def));

fidR = fopen("../../data/holChirpRealMatlab.dat",'w');
fwrite(fidR,phase_diff,'double');
fclose(fidR);

break;

%MS_Real = real(chirp);
%MS_Imag = imag(chirp);
%fidR = fopen("../../data/chirpReal.dat",'w');
%fwrite(fidR,MS_Real,'double');
%fclose(fidR);
%fidI = fopen("../../data/chirpImag.dat",'w');
%fwrite(fidI,MS_Imag,'double');
%fclose(fidI);

status = sprintf('done with chirp. \n')

status = sprintf('Computing Fresnel transform of def...')
b_def =fft2(hol_def.*chirp);
b_def =fftshift(b_def);
%b_def = ifftshift(ifft2(fftshift(hol_def.*chirp)));
%status = sprintf('done with Fresnel transform of def. \n')

%status = sprintf('Computing Fresnel transform of base...')
b_base =fft2(hol_base.*chirp);
b_base =fftshift(b_base);
%b_base = ifftshift(ifft2(fftshift(hol_base.*chirp)));

fidBR = fopen("../../data/baseReal.dat",'w');
fwrite(fidBR,real(b_base),'double');
fclose(fidBR);
fidBI = fopen("../../data/baseImag.dat",'w');
fwrite(fidBI,imag(b_base),'double');
fclose(fidBI);

status = sprintf('done with Fresnel transform of base. \n')


status = sprintf('Computing interference phase...')
phase_diff = atan2(imag(b_base),real(b_base));
status = sprintf('done computing interference phase. \n')

fnt = 20;
xsize_fres = abs(d*lambda/M/res);
ysize_fres = abs(d*lambda/N/res);
x = 0:xsize_fres:(M-1)*xsize_fres;
y = 0:ysize_fres:(N-1)*ysize_fres;
x_plt = x - max(x)/2;
y_plt = y - max(y)/2;
    
figure(1); hold on;
fig1 = imagesc(x,y,phase_diff); axis image;%colormap gray; 
set(gca,'Ydir','normal');
title('Reconstructed Phase [rad]','fontsize',fnt);
ylabel('\xi [m]','fontsize',fnt);
xlabel('\eta [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
colorbar;




chirp = ones(N,M);

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

L = size(hol_base);
M = L(2);  % x-direction
N = L(1);  % y-direction

param = 8;
x2 = linspace(-param*pi,param*pi,M+1); x = x2(1:M);
y2 = linspace(-param*pi,param*pi,N+1); y = y2(1:N);

for i = 1: length(hol_base(:,1))
    hol_def(i,:) = hol_def(i,:).*(tanh(x+(param-1)*pi)-tanh(x-(param-1)*pi));
    hol_base(i,:) = hol_base(i,:).*(tanh(x+(param-1)*pi)-tanh(x-(param-1)*pi));
end

for j = 1: length(hol_base(1,:))
    hol_def(:,j) = hol_def(:,j).*(tanh(y+(param-1)*pi)-tanh(y-(param-1)*pi))';
    hol_base(:,j) = hol_base(:,j).*(tanh(y+(param-1)*pi)-tanh(y-(param-1)*pi))';
end	



end
