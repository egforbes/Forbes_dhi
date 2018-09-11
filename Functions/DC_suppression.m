% DC_suppression.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code suppresses the zeroth order diffraction image by
% subtracting the mean intensity from each hologram.


function [hol_def,hol_base] = DC_suppression(hol_def,hol_base,M,N)

hol_def = hol_def - ones(N,M)*mean(mean(hol_def));
hol_base = hol_base - ones(N,M)*mean(mean(hol_base));

status = sprintf('DC diffraction image suppressed... \n')