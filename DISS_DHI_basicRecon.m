% DHI_basicRecon.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs the Fresnel transform reconstruction to
% obtain the phase difference between specified plasma and vacuum
% holograms.  This code primarily allows the user to view the
% reconstruction plane to identify the bounds of the twin images for
% further analysis.  The code proceeds through 9 main sections:
%
% 1.) Declare reconstruction parameters/constants
% 2.) Input holograms from file
% 3.) Mask filter the hologram
% 4.) DC image suppression
% 5.) Apply windowing to the holograms
% 6.) Reference wave generation
% 7.) Fresnel transform reconstruction
% 8.) Compute the phase difference
% 9.) Plotting the reconstruction plane

tic; clear all; clc;
close all; drawnow;
format compact;

section = 0;
cross_sect = section;

shotnum_base_vec = [180621015]; % user specifies vacuum/baseline shot here

for shotnum_def_it = [180621015]; % user specifies plasma shot here
    close all, drawnow;
    close all, drawnow;
    
    ind_base_tmp = find((shotnum_def_it-shotnum_base_vec)>0,1,'last');
%     shotnum_base_it = shotnum_base_vec(ind_base_tmp);
shotnum_base_it = shotnum_base_vec(1);
    
    %% 1.) Declare reconstruction parameters:
    % Significantly, the reconstruction
    % distance, d, only affects the resolution of the reconstructed
    % interference phase.  Maximize resolution by selecting the 
    % reconstruction distance to maximize the size of the twin images while
    % still separating them from the DC term.
    
    status = sprintf('Initializing... \n')
    
    res = (3.85e-6);% [m] Pixel size of Nikon 3200D, 4.3e-6 for Canon T2i
    lambda = 532e-9;% [m] Nd:YAG wavelength, ruby is 694.3e-9 [m]
    d =.4; % [m] Reconstruction distance
    n_c = 2.312e27; % [m^-3] Critical plasma density
    
    %% 2.) Input holograms from file:
    % This section uses imread() to load MxN-sized uint8 datatype color
    % holograms into MATLAB.  Each MxN color hologram is stored in an MxNx3
    % array.  Each of the three MxN slices correspond to the intensity
    % distribution of each color in the the RGB color palette.  To consider
    % the total intensity at each pixel, this code converts the MxNx3 array
    % to an MxN grayscale array with rgb2gray().  The grayscale intensity 
    % value at a given pixel location amounts to a norm of the three RGB 
    % intensities at that pixel.  Finally, the code converts the uint8 
    % datatype to double to prepare the intensity information for 
    % mathematical manipulation.
    
    %%% Starting with ZaP's shot number convention:
    SN_def = num2str(shotnum_def_it);
    SN_base = num2str(shotnum_base_it);
    
    %%% Converting to Canon Camera shot number convention:
    year_def = num2str(2000+str2num(SN_def(1:2)));
    month_def = num2str(SN_def(3:4));
    day_def = num2str(SN_def(5:6));
    
    year_base = num2str(2000+str2num(SN_base(1:2)));
    month_base = num2str(SN_base(3:4));
    day_base = num2str(SN_base(5:6));
    
    date_def = [year_def,'_', month_def,'_', day_def];
    date_base = [year_base,'_',month_base,'_',day_base];
    date_plot = [year_def, month_def, day_def];
    
    shot_def = SN_def(7:end);
    shot_base = SN_base(7:end);
    
    % Specifying where the raw holograms are stored:
%      hol_def_rgb = imread(['M:\Users\Eleanor_Forbes\',date_def,'\',...
%         date_def, '_', shot_def,'.JPG']);
%     hol_base_rgb = imread(['M:\Users\Eleanor_Forbes\',date_base,'\',...
%         date_base, '_', shot_base,'.JPG']);
hol_def_rgb = imread(['K:\Fuze Data\DHI_Images\','180621','\',...
    '180621', shot_def,'.JPG']);
hol_base_rgb = imread(['K:\Fuze Data\DHI_Images\','180621','\',...
    '180621', shot_def,'_Baseline.JPG']);
    
    
%     hol_def_rgb = imread(['H:\Raw Holograms\',date_def,'\',...
%         date_def, '_', shot_def,'.jpg']);
%     hol_base_rgb = imread(['H:\Raw Holograms\',date_base,'\',...
%         date_base, '_', shot_base,'.jpg']);
%     
    % Converting from RGB to grayscale (rgb2gray is not a supported command
    % without the image processing toolbox):
    rgb_red = 0.3;
    rgb_green = 0.59;
    rgb_blue = 0.11;
    
    hol_def = double(hol_def_rgb(:,:,1).*rgb_red + hol_def_rgb(:,:,2).*rgb_green + hol_def_rgb(:,:,3).*rgb_blue)';
    hol_base = double(hol_base_rgb(:,:,1).*rgb_red + hol_base_rgb(:,:,2).*rgb_green + hol_base_rgb(:,:,3).*rgb_blue)';
    
    %% 3.) Mask filter the hologram:
    % This section applies a mask filter by selecting only a limited 
    % section of the grayscale holograms.  The mask filter can improve the 
    % quality of the reconstruction by cutting out unused parts of the 
    % camera sensor.
    
    % xmin = 1200;
    % xmax = 4200;
    % ymin = 1000;
    % ymax = 4000;
    % hol_def = hol_def(ymin:ymax,xmin:xmax);
    % hol_base = hol_base(ymin:ymax,xmin:xmax);
    
    % Determining the size of the holograms:
    L = size(hol_def);
    M = L(2);  % x-direction
    N = L(1);  % y-direction
    
    %% 4.) DC image suppression:
    % This section subtracts the mean intensity of the hologram from the
    % intensity at each pixel location.  In some cases, this can improve 
    % the reconstructed intensity.  For the purpose of reconstructing the 
    % phase, the benefit is minimal as long as the twin images are 
    % sufficiently separated from the DC term, so this section can be 
    % omitted to save computation time. Subtracting the mean intensity from
    % the hologram serves to suppress the DC diffraction term by 
    % eliminating the DC (or zeroth) order term in the discrete Fourier 
    % series of the hologram. Our reconstruction relies on the discrete 
    % Fourier transform, so removing its DC term reduces the DC diffraction
    % term, improving the reconstructed intensity.  Recall that the DC term
    % of the Fourier series is merely the mean value of the function.
    
    % [hol_def,hol_base] = DC_suppression(hol_def,hol_base,M,N);
    
    %% 5.) Apply windowing to the hologram:
    % Applying a window to the raw holograms before reconstruction 
    % can reduce Gibb's phenomena by smoothing out non-periodic boundary
    % conditions.  This provides marginal noise reduction in the 
    % reconstructed interference phase distribution.  Two window functions 
    % are available here.
    
    % Hamming window:
    % [hol_def,hol_base] = hamming_window(hol_def,hol_base);
    % Hyperbolic window:
    % [hol_def,hol_base] = hyperbolic_window(hol_def,hol_base);
    
    %% 6.) Reference wave generation:
    % Because our hologram recording utilizes planar reference waves, we 
    % can model our reference beam as a plane wave.  If we use lenses to 
    % create a spherical reference beam during recording, modeling the 
    % reference wave will involve a complex exponential.
    
    ref = 1;
    
    %% 7.) Fresnel transform reconstruction:
    % The Fresnel transform reconstruction method generates the electric
    % field strength at reconstruction distance d for the deformed and 
    % baseline holograms.  The Fresnel method, which assumes parabolic 
    % wavelets as source functions, is an approximation of scalar 
    % diffraction theory applicable for sufficiently large reconstruction 
    % distances.
    
    b_def = fresnel(hol_def,d,lambda,res,ref);
    b_base = fresnel(hol_base,d,lambda,res,ref);
    
    %% 8.) Compute the phase difference:
    % Properties of complex numbers are applied to the complex wavefield to
    % compute the phase difference.
    
    phase_diff = phase_diff_func(b_def,b_base);
    
    %% 9.) Plotting the reconstruction plane:
    
    fnt = 20;
    xsize_fres = abs(d*lambda/M/res);
    ysize_fres = abs(d*lambda/N/res);
    x = 0:xsize_fres:(M-1)*xsize_fres;
    y = 0:ysize_fres:(N-1)*ysize_fres;
    x_plt = x - max(x)/2;
    y_plt = y - max(y)/2;
    
    figure(1); hold on;
%     fig1 = imagesc(x,y, phase_diff); axis image; %colormap gray; 
    fig1 = imagesc(flipud(fliplr(phase_diff.'))); axis image
    set(gca,'Ydir','normal');
    title('Reconstructed Phase [rad]','fontsize',fnt);
    ylabel('\xi [m]','fontsize',fnt);
    xlabel('\eta [m]','fontsize',fnt);
    set(gca,'fontsize',fnt);
    colorbar;
    
end