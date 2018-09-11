% phase_diff_func.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code computes the wrapped phase difference between
% plasma and vacuum holograms.

function [phase_diff] = phase_diff_func(b_def,b_base)

status = sprintf('Computing interference phase...')
phase_diff = atan2(imag((b_def).*conj((b_base))),real((b_def).*...
    conj((b_base))));
status = sprintf('done computing interference phase. \n')