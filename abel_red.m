% abel_invert.m

function [den_int,Anew] = abel_red(radius,den_num,abel_diff)

dr = radius(2)-radius(1);

sizeold = length(den_num)+abel_diff;

A = zeros(sizeold,sizeold);

for k = sizeold:-1:1
    for i = sizeold:-1:1
       
        A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
        
    end
end


% A = 2*dr*real(A');
% 
% den_int = A*den_num';


A = A';
sizenew = length(den_num);
Anew = zeros(sizenew,sizenew);

for kk = 1:sizenew

    Anew(kk,:) = A(kk,1:sizenew);
    Anew(kk,sizenew) = sum(A(kk,sizenew:end));

end

Anew = 2*dr*real(Anew);

den_int = Anew*den_num';





% function [den_int] = abel(radius,den_num)
% 
% dr = radius(2)-radius(1);
% 
% A = zeros(length(den_num),length(den_num));
% 
% for k = length(den_num):-1:1
%     for i = length(den_num):-1:1
%        
%         A(k,i) = sqrt(((k+1))^2-(i)^2)-sqrt((k)^2-(i)^2);
%         
%     end
% end
% 
% A = 2*dr*real(A');
% 
% den_int = A*den_num';
