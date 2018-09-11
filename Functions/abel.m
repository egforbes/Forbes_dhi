% abel_invert.m

function [den_int] = abel_invert(radius,den_num)

dr = radius(2)-radius(1);

A = zeros(length(den_num),length(den_num));

for k = length(den_num):-1:1
    for i = k:-1:1
       
        A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
        
    end
end

A = 2*dr*(A');

den_int = A*den_num';

