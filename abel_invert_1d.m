% abel_invert.m

function [den_num,Anew] = abel_invert_1d(radius,den_int,abel_diff)

dr = radius(2)-radius(1);

sizeold = length(den_int)+abel_diff;

A = zeros(sizeold,sizeold);

for k = sizeold:-1:1
    for i = sizeold:-1:1
       
        A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
        
    end
end

A = A';
sizenew = length(den_int);
Anew = zeros(sizenew,sizenew);

for kk = 1:sizenew

    Anew(kk,:) = A(kk,1:sizenew);
    Anew(kk,sizenew) = sum(A(kk,sizenew:end));

end

Anew = 2*dr*real(Anew);

den_num = Anew\den_int';







% function [den_num,A] = abel_invert_1d(radius,den_int)
% 
% dr = radius(2)-radius(1);
% 
% A = zeros(length(den_int),length(den_int));
% 
% for k = length(den_int):-1:1
%     for i = length(den_int):-1:1
%        
%         A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
%         
%     end
% end
% 
% A = 2*dr*real(A');
% 
% den_num = A\den_int';
