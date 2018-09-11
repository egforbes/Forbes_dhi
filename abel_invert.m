% abel_invert.m

function [den_num] = abel_invert(radius,den_int)

dr = radius(2)-radius(1);

A = zeros(length(den_int),length(den_int));

for k = length(den_int):-1:1
    for i = k:-1:1
       
        A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
        
    end
end

A = 2*dr*(A');

den_num = A\den_int;
