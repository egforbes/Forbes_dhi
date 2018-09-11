% holographic_reconstruction.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs the Fresnel transform reconstruction and
% computes the phase shift between specified plasma and vacuum holograms.
% It also extracts a twin image from the data and smooths and unwraps that
% twin image.

function [unwrapped,x_twin,y_twin,phase_diff,x,y] = ...
    holographic_reconstruction(shotnum_base,shotnum_def,xmin,xmax,ymin,...
    ymax,d,varargin)

% Constants:
res = (3.85e-6); % [m] Pixel size of Nikon 3200D 4.3e-6 for Canon Rebel T2i

% Starting with ZaP's shot number convention:
SN_def = num2str(shotnum_def);
SN_base = num2str(shotnum_base);

% Converting to Canon Camera shot number convention:
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

hol_def_rgb = imread(['M:\M_drive\Users\Eleanor_Forbes\',date_def,'\',...
    date_def, '_', shot_def,'.JPG']);
hol_base_rgb = imread(['M:\M_drive\Users\Eleanor_Forbes\',date_base,'\',...
    date_base, '_', shot_base,'.JPG']);
% 
% hol_def_rgb = imread(['M:\Users\Eleanor_Forbes\',date_def,'\',...
%     date_def, '_', shot_def,'.JPG']);
% hol_base_rgb = imread(['M:\Users\Eleanor_Forbes\',date_base,'\',...
%     date_base, '_', shot_base,'.JPG']);
% hol_def_rgb = imread(['K:\Fuze Data\DHI_Images\','180710','\',...
%     '18071', shot_def,'.JPG']);
% hol_base_rgb = imread(['K:\Fuze Data\DHI_Images\','180710','\',...
%     '18071', shot_def,'_Baseline.JPG']);
% hol_def_rgb = imread(['K:\Fuze Data\DHI_Images\Calibration\',...
%     'June_25_2018\Flame_2.JPG']);
% hol_base_rgb = imread(['K:\Fuze Data\DHI_Images\Calibration\',...
%     'June_25_2018\NoFlame_1.JPG']);
%     

% hol_def_rgb = imread(['H:\Raw Holograms\',date_def,'\',...
%     date_def, '_', shot_def,'.jpg']);
% hol_base_rgb = imread(['H:\Raw Holograms\',date_base,'\',...
%     date_base, '_', shot_base,'.jpg']);

% Converting from RGB to grayscale (rgb2gray is not a supported command
% without the image processing toolbox):
rgb_red = 0.3;
rgb_green = 0.59;
rgb_blue = 0.11;

hol_def = double(hol_def_rgb(:,:,1).*rgb_red + hol_def_rgb(:,:,2).*...
    rgb_green + hol_def_rgb(:,:,3).*rgb_blue)';
hol_base = double(hol_base_rgb(:,:,1).*rgb_red + hol_base_rgb(:,:,2).*...
    rgb_green + hol_base_rgb(:,:,3).*rgb_blue)';

% If an eighth input = 1 is passed to the function, the hologram intensity
% values are truncated within the range hol_max_cut:hol_min_cut.
if nargin == 8
    hol_contrast = varargin{1};
    
    hol_max_cut = 100;
    hol_min_cut = 70;
    
    if hol_contrast == 1;
        hol_def(find(hol_def>hol_max_cut)) = hol_max_cut;
        hol_def(find(hol_def<hol_min_cut)) = hol_min_cut;
    else
    end
else
end

% Reconstruction:
[phase_diff,x,y,xsize_fres,ysize_fres] = DHI_rec(hol_base,hol_def,d,res);

M = size(hol_def,2);
N = size(hol_def,1);

x_plt = -max(x)/2:x(2)-x(1):max(x)/2;
y_plt = -max(y)/2:y(2)-y(1):max(y)/2;

% Extracting the twin image:

% Labeling the matrix indices:
x_ind_vec = 1:1:M;
y_ind_vec = 1:1:N;
 
% Interpolating to match a spatial coordinate values to the
% nearest integer matrix indices:
x_min_ind = round(interp1(x,x_ind_vec,xmin));
x_max_ind = round(interp1(x,x_ind_vec,xmax));
y_min_ind = round(interp1(y,y_ind_vec,ymin));
y_max_ind = round(interp1(y,y_ind_vec,ymax));

% Defining the x,y spatial coordinate vectors for the extracted twin image:
x_twin = x(x_min_ind:x_max_ind);
x_twin = x_twin(1:2:end);
y_twin = y(y_min_ind:y_max_ind);
y_twin = y_twin(1:2:end);
% 
twin_img = phase_diff(y_min_ind:y_max_ind,x_min_ind:x_max_ind);
twin_img = twin_img(1:2:end,1:2:end);

% Smoothing and unwrapping the twin image:
% width = 80;
width = 30; %8 for FuZE, 4 for HD
[unwrapped,phase_f,unsmoothed_unwrapped] = smooth_unwrap(twin_img,width);