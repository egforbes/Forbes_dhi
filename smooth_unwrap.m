% smooth_unwrap.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code smooths the extracted wrapped twin image and then
% unwraps the phase.  The function outputs smoothed and unsmoothed versions
% of the unwrapped phase.

function [unwrapped,phase_f,unsmoothed_unwrapped] = smooth_unwrap(...
    twin_img,width)

% Transforming into cosine/sine space:
s = sin(twin_img); 
c = cos(twin_img); 

% Boxcar averaging (change the width parameter to change the size of the
% sliding boxcar window):
% width = 4;
sf = boxcar2(s,width); 
cf = boxcar2(c,width); 

% Inverting back to "phase space":
phase_F = atan2(sf,cf);
phase_f = phase_F(2*width+1:end-(2*width+1),:);

% unwrapping the smoothed extracted twin image:
unwrapped = zeros(size(phase_f,1),size(phase_f,2));

% unwrap rows
for i =1:size(phase_f,1)
    unwrapped(i,:) = unwrap(phase_f(i,:));
end
% unwrap columns
for j =1:size(phase_f,2) 
    unwrapped(:,j) = unwrap(unwrapped(:,j));
end

unwrapped = unwrapped - max(unwrapped(:));  

% unwrapping the unsmoothed extracted twin image:
unsmoothed_unwrapped = zeros(size(twin_img,1),size(twin_img,2));

% unwrap rows
for i =1:size(phase_f,1)
    unsmoothed_unwrapped(i,:) = unwrap(twin_img(i,:));
end
% unwrap columns
for j =1:size(phase_f,2) 
    unsmoothed_unwrapped(:,j) = unwrap(unsmoothed_unwrapped(:,j));
end

unsmoothed_unwrapped = unsmoothed_unwrapped - max(unsmoothed_unwrapped(:));  