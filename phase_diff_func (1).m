%% Compute phase difference between baseline and plasma shot:

function [phase_diff] = phase_diff_func(b_def,b_base)

status = sprintf('Computing interference phase...')
phase_diff = atan2(imag((b_def).*conj((b_base))),real((b_def).*conj((b_base))));
status = sprintf('done computing interference phase. \n')