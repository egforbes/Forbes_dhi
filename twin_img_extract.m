function [x_twin,y_twin,twin_img] = twin_img_extract(x,y,phase_diff,x_hol,y_hol,h)

L = size(phase_diff);
M = L(2);  % x-direction
N = L(1);  % y-direction

%% Extracting the twin image:

% Labeling the matrix indices:
x_ind_vec = 1:1:M;
y_ind_vec = 1:1:N;

% Identifying the confines of the twin image:
xmin = -max(x_hol)/2;
xmax = max(x_hol)/2;
ymin = h - max(y_hol)/2;
ymax = h + max(y_hol)/2;
% testing sensitivity to changing mask size:
% xmin = -max(x_hol-1e-5)/2;
% xmax = max(x_hol-1e-5)/2;
% ymin = h - max(y_hol-1e-5)/2;
% ymax = h + max(y_hol-1e-5)/2;

% Interpolating to match a spatial coordinate values to the
% nearest integer matrix indices:
x_min_ind = round(interp1(x,x_ind_vec,xmin));
x_max_ind = round(interp1(x,x_ind_vec,xmax));
y_min_ind = round(interp1(y,y_ind_vec,ymin));
y_max_ind = round(interp1(y,y_ind_vec,ymax));

% Defining the x,y spatial coordinate vectors for the extracted twin image:
x_twin = x(x_min_ind:x_max_ind);
y_twin = y(y_min_ind:y_max_ind);

% Extracting the twin image:
twin_img = phase_diff(y_min_ind:y_max_ind,x_min_ind:x_max_ind);

% Centering the spatial coordinates for the twin image:
x_twin = x_twin - x_twin(round(length(x_twin)/2));
y_twin = y_twin - y_twin(round(length(y_twin)/2));
