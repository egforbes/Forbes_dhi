clear all; clc;
close all;

%% Ross's uniform drift vel:
clear all; clc;
close all;
k = 1;
mu0=4e-7*pi; %N/A^2
shotnum = 160524021;
datafile = ['dhi_error_data.mat'];
filepath='M:\Users\Michael_Ross\Coding\Holography_Reconstruction\';
filename=[filepath,datafile]
load(filename);
r=rad_l{k}';
n=den_num_axial_l_plt{k};

% r = linspace(0,0.02,100);
% a = 0.01;
% n_0_2 = 5e22;
% n = [n_0_2*(a^2./(r.^2+a^2))];
% I = 40e3;

[time,m_0_p0] = MDS_MATLAB_python(shotnum,'m_0_p0');
[time_holo] = MDS_MATLAB_python_digitizers(shotnum,...
    'ZAPHD::TOP.DIGITIZERSHD.L4222:ZAP_4222_B:CHANNEL_1');
ind_time = find((time_holo - time)<0,1,'first');
Bwall = m_0_p0(ind_time);
rwall = 4*2.54/100;
Iz = Bwall*2*pi*rwall/mu0;

[B_prof,T_prof,v_drift,den_char,B_char,T_char,a_char] = DHI_profiles(r,n,Iz,.01);

figure(1);
subplot(3,1,1); hold on;
plot(r,n);
subplot(3,1,2); hold on;
plot(r,B_prof);
subplot(3,1,3); hold on;
plot(r,T_prof);

k_boltz = 1.38e-23 %J/K; %8.62e-5;% eV/K
T_prof = 11600*T_prof;

n_red = n(1:length(n)-1);
T_red = T_prof(1:length(T_prof)-1);
B_red = B_prof(1:length(B_prof)-1);
r_red = r(1:length(r)-1);

LHS = (2*k_boltz*n_red.*(diff(T_prof)./diff(r)) + 2*k_boltz*T_red.*(diff(n)./diff(r)));
RHS = (-B_red.^2./(mu0*r_red) - (B_red/mu0).*(diff(B_prof)./diff(r)));
% RHS2 = (-B_red./(mu0*r_red)).*(diff(r.*B_prof)./diff(r));

% figure(2);
% subplot(4,1,1);
% plot(diff(T_prof));
% subplot(4,1,2);
% plot(diff(n));
% subplot(4,1,3);
% plot(diff(B_prof));
% subplot(4,1,4);
% plot(diff(r));

figure(3);
plot(LHS./RHS);

%% Stu's uniform drift vel:
clear all; clc;
% shotnum = 160524021;
% % datafile = ['dhi_error_data_',num2str(shotnum),'.mat'];
% datafile = ['dhi_error_data.mat'];
% a = stuart_equil(datafile,shotnum)

k = 1;
mu0=4e-7*pi; %N/A^2
shotnum = 160524021;
datafile = ['dhi_error_data.mat'];
filepath='M:\Users\Michael_Ross\Coding\Holography_Reconstruction\';
filename=[filepath,datafile]
load(filename);
r_full=rad_l{k}';
n_full=den_num_axial_l_plt{k};

[B,T,v_drift,den_char,B_char,T_char,a_char] = stuart_equil(r_full,n_full,shotnum);

figure(101);
subplot(3,1,1); hold on;
plot(r_full,n_full);
scatter(r_full(1),den_char);
subplot(3,1,2); hold on;
plot(B);
scatter(length(B),B_char);
subplot(3,1,3); hold on;
plot(T);
scatter(1,T_char);

k_boltz = 1.38e-23;
T = T*11600;

n_red = n_full(1:length(n_full)-1);
T_red = T(1:length(T)-1);
B_red = B(1:length(B)-1);
r_red = r_full(1:length(r_full)-1);

LHS = 2*k_boltz*n_red.*(diff(T)./diff(r_full)) + 2*k_boltz*T_red.*(diff(n_full)./diff(r_full));
RHS = -B_red.^2./(mu0*r_red) - (B_red/mu0).*(diff(B)./diff(r_full));
% RHS2 = (-B_red./(mu0*r_red)).*(diff(r_full.*B)./diff(r_full));

% figure(102);
% subplot(4,1,1);
% plot(diff(T));
% subplot(4,1,2);
% plot(diff(n_full));
% subplot(4,1,3);
% plot(diff(B));
% subplot(4,1,4);
% plot(diff(r_full));

figure(103);
plot(LHS./RHS);

%% Stu's uniform temperature:
% shotnum = 160524021;
% datafile = ['dhi_error_data.mat'];
% % [B,T,den_char,B_char,T_char,a_char] = temp_calc_const_t(datafile,shotnum);
% 
% 
% k = 1;
% filepath='M:\Users\Michael_Ross\Coding\Holography_Reconstruction\';
% filename=[filepath,datafile]
% load(filename);
% 
% r_full = rad_l{k}';
% n_full=den_num_axial_l_plt{k};
% [B,T,v_drift,den_char,B_char,T_char,a_char]=temp_calc_const_t(r_full,n_full,shotnum)

