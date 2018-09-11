
% error_in is N_error_l{k}
% radius is rad_rec_l{k}

function [n_error_out] = abel_invert_uncertainty(radius,N_error_in)


dr = radius(2)-radius(1);

A = zeros(length(N_error_in),length(N_error_in));

for k = length(N_error_in):-1:1
    for i = k:-1:1
       
        A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
        
    end
end

A = (A');

n_error_out = ((1/(2*dr))^2*inv(A).^2*N_error_in'.^2).^(0.5);


