% boxcar2.m

% Description:  This code conducts a two-dimensional boxcar smoothing on the
% two-dimensional dataset contained in matrix "A."  The algorithm iterates a
% Shannon (square) window across the two-dimensional data.  At each
% iteration, the window (of side length 2*width+1) is centered at a
% different element in matrix "A."  The value of this element is replaced by
% the mean value of all the elements contained window.  Note that the size
% of the Shannon window can be modified by changing the "width" parameter.

% Inputs:
% A:  Two-dimensional dataset
% width:  scalar value used to set the width of the Shannon window

% Outputs:
% smoothed:  Two-dimensional boxcar-smoothed dataset

function [smoothed] = boxcar2(A,width)

l = size(A);
m = l(2); % x-direction
n = l(1); % y-direction

for i = 1:m*n

    [I,J] = ind2sub(l,i);
    
    I_tmp = [I-width:I+width];
    J_tmp = [J-width:J+width];
    
    % At the edges of the dataset, we forgo averaging and simply keep the
    % unsmoothed values of the dataset. 
    I_ind = find((I_tmp > 0) & (I_tmp < size(A,1)));
    J_ind = find((J_tmp > 0) & (J_tmp < size(A,2)));
    
    if length(I_ind)<(2*width+1)
        smoothed(I,J) = A(I,J);
    else
        avg_tmp = A(I_tmp(I_ind),J_tmp(J_ind));
        smoothed(I,J) = mean(avg_tmp(:));
    end
    
end
