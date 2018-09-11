% DHI_rec.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs the Fresnel transform reconstruction and
% computes the phase shift between specified plasma and vacuum holograms.

function [phase_diff,x,y,xsize_fres,ysize_fres] = DHI_rec(hol_base,...
    hol_def,d,res)

%% Declare parameters:
lambda = 532e-9; % [m] Nd:YAG wavelength, 694.3e-9 [m] ruby wavelength 
n_c = 2.312e27; % [m^-3] critical plasma density

L = size(hol_def);
M = L(2);  % x-direction
N = L(1);  % y-direction

ximg = linspace(0,res*M,M);
yimg = linspace(0,res*N,N);

%% Applying windowing to the holograms:
[hol_def,hol_base] = hyperbolic_window(hol_def,hol_base);
% [hol_def,hol_base] = DC_suppression(hol_def,hol_base,M,N);

%% Reconstruction with Fresnel method:
ref = 1;
b_def = fresnel(hol_def,d,lambda,res,ref);
b_base = fresnel(hol_base,d,lambda,res,ref);

%% Computing the phase from the complex wavefield:
phase_diff = phase_diff_func(b_def,b_base);

%% Preparing axes necessary for plotting:
xsize_fres = abs(d*lambda/M/res);
ysize_fres = abs(d*lambda/N/res);
x = 0:xsize_fres:(M-1)*xsize_fres;
y = 0:ysize_fres:(N-1)*ysize_fres;