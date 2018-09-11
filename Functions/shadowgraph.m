function [dn,dN,A_test,del_I,del_I_mat] = shadowgraph(radius,n_profile,L,M,N)





lambda = 533e-9;
q_e = 1.6e-19;
c = 3.0e8;
m_e = 9.12e-31;
eps_0 = 8.85e-12;

dn = 0.5*n_profile*lambda^2*q_e^2/(c^2*m_e*eps_0*4*pi^2);

dN = abel(radius,dn);


n = length(radius);

X_0 = -2*ones(n,1);
X_1 = ones(n,1);

LHB = zeros(n,1);
RHB = zeros(n,1);
LHB(1:4) = [2 -5 4 -1];
RHB(n-3:n) = [-1 4 -5 2];

Atest = spdiags([X_1 X_0 X_1],[-1 0 1],n,n);
Atest(1,:) = LHB;
Atest(end,:) = RHB;

dR = radius(2)-radius(1);
A = spdiags([X_1 X_0 X_1],[-1 0 1],n,n)*(1/(dR^2));
A(1,:) = LHB*(1/(dR^3));
A(end,:) = RHB*(1/(dR^3));

A_test = full(A);
% figure(2);
% spy(A);


%%% Solving for a relative intensity distribution:
del_I = L*A*dN;
del_I(1) = del_I(2);
del_I(end) = del_I(length(del_I)-1);

radius_lg = linspace(0,max(radius),M);
del_I_mat_seed = interp1(radius,del_I,radius_lg);
del_I_mat = ones(M,N);
for i = 1:M
    del_I_mat(i,:) = del_I_mat_seed;
end
