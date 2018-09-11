function [rad_rec,den_int_rec,den_num_rec,ph_rec,A] = phase2profiles(unwrapped,xsize_fres,ysize_fres,abel_diff,N_edge,l)
% function [rad_rec,rad_s_rec,den_int_rec,den_num_rec,ph_rec,den_int_rec_s] = phase2profiles(unwrapped,xsize_fres,ysize_fres)

% Defining constants:
e_0 = 8.85e-12;
mu_0 = 1.257e-6;
k = 1.38e-23; % Boltzmann's constant in SI units
q = 1.6e-19; % electron charge [C]
c = 3e8; % speed of light [m/s]
mi = 1.67e-27; % proton mass [kg]
me = 9.12e-31; % electron mass [kg]
lambda = 532e-9;

% Orienting the reconstructed phase to match the orientation of the assumed
% density profile:
ph_rec = -fliplr(unwrapped(round(size(unwrapped,1)/2),:)); % note sign here depends on which first order term used

% Converting phase to line-integrated density:
den_int_rec = ph_rec/((-q^2/(4*pi*c^2*me*e_0))*lambda);
den_int_rec = (den_int_rec-min(den_int_rec));%   for the adjusted Abel inversion method:  +N_edge;

% Defining radial dimension vectors:
rad_rec = [0:xsize_fres:((size(unwrapped,2)-1)*xsize_fres)];   
% rad_s_rec= [0:200e-6:((size(unwrapped,2)-1)*xsize_fres)];

% Resampling the line-integrated density onto a sparser grid (this saves
% tons of time during the Abel inversion:
% den_int_rec_s = interp1(rad_rec,den_int_rec,rad_s_rec);

% Obtaining number density through Abel inversion of the line-integrated
% density:
% den_num_rec = abel_invert_1d(rad_s_rec,den_int_rec_s);
if l ==0
    [den_num_rec,A] = abel_invert_1d(rad_rec,den_int_rec,abel_diff);
else
    [den_num_rec,A] = abel_invert_1d(rad_rec,den_int_rec,abel_diff);
    den_num_rec = den_num_rec + N_edge/l;
end