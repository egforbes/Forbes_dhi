% hyperbolic_window.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs applies a hyperbolic window function to
% the holograms to remove Gibbs phenomena in the resulting recontructions.

function [hol_def, hol_base] = hyperbolic_window(hol_def,hol_base)

L = size(hol_def);
M = L(2);  % x-direction
N = L(1);  % y-direction

param = 8;
x2 = linspace(-param*pi,param*pi,M+1); x = x2(1:M);
y2 = linspace(-param*pi,param*pi,N+1); y = y2(1:N);

for i = 1: length(hol_def(:,1))
    hol_def(i,:) = hol_def(i,:).*(tanh(x+(param-1)*pi)-tanh(x-(param-1)*pi));
    hol_base(i,:) = hol_base(i,:).*(tanh(x+(param-1)*pi)-tanh(x-(param-1)*pi));
end

for j = 1: length(hol_base(1,:))
    hol_def(:,j) = hol_def(:,j).*(tanh(y+(param-1)*pi)-tanh(y-(param-1)*pi))';
    hol_base(:,j) = hol_base(:,j).*(tanh(y+(param-1)*pi)-tanh(y-(param-1)*pi))';
end