% clear all; close all; clc;

close all

%PRL Plot generation

load prl_dhi_err3.mat;
%structures from this file are 'rad_l','n_prof_l','rad_r','n_prof_r','n_error_l','n_error_r'
indd = 1;
figure(1)
plot(rad_l{indd},n_prof_l{indd})
hold on
plot(rad_r{indd},n_prof_r{indd})

if length(rad_l{indd}(1:end)) < length(rad_r{indd}(1:end))
    r_ax = rad_l{indd}(1:end);
    l_prof = n_prof_l{indd};
    r_prof = n_prof_r{indd}(1:length(l_prof));
    l_err = n_error_l{indd};
    r_err = n_error_r{indd}(1:length(l_prof));
%     B_l = B_prof_l{cross_sect(indd)};
%     B_r = B_prof_r{cross_sect(indd)}(1:length(l_prof));
%     T_l = T_prof_l{cross_sect(indd)};
%     T_r = T_prof_r{cross_sect(indd)}(1:length(l_prof));
%     
else
    r_ax = rad_r{indd}(1:end);
    r_prof = n_prof_r{indd};
    l_prof = n_prof_l{indd}(1:length(r_prof));
    l_err = n_error_l{indd}(1:length(r_prof));
    r_err = n_error_r{indd};
%     B_r = B_prof_r{cross_sect(indd)};
%     B_l = B_prof_l{cross_sect(indd)}(1:length(r_prof));
%     T_r = T_prof_r{cross_sect(indd)};
%     T_l = T_prof_l{cross_sect(indd)}(1:length(r_prof));
    
end

n_mean = (l_prof+r_prof)./2E6; %take average and convert to cm^-3
err_mean = (abs(l_err(1:end-1))+abs(r_err(1:end-1)))./2E6; %take average and convert to cm^-3
% T_mean = (T_r + T_l)./2;
% B_mean = (B_r+B_l)./2;

figure(2)
plot(r_ax,n_mean,'-k','LineWidth',2)
hold on
errorbar(r_ax(1:8:end-1),n_mean(1:8:end-1),err_mean(1:8:end),'ko','LineWidth',2)
set(gca,'Fontsize',14)
xlabel('Impact Parameter (cm)')
set(gca,'xtick',[0:2E-3:6E-3])
set(gca,'xticklabel',{'0','0.2','0.4','0.6'})
xlim([-0.1E-3 6E-3])
ylabel('n_e (cm^{-3})')
title('Number Density at z = 10 cm')

mdsconnect('zappa.zap');
mdsopen('zaphd',171005023);
m_0_p0 = mdsvalue('\m_0_p0');
time = mdsvalue('dim_of(\m_0_p0)');
time_holo = 45E-6;

mu0 = (4*pi)*1e-7;
ind_time = find((time_holo - time)<0,1,'first');
Bwall = m_0_p0(ind_time);
rwall = 4*2.54/100;
Iz = Bwall*2*pi*rwall/mu0;
mdsclose
mdsdisconnect

% [B_prof,T_prof,v_drift,den_char,B_char,T_char,a_char] = DHI_profiles(r_ax,n_mean,Iz);
% [B_err,T_err,v_err,den_c_err,B_c_err,T_c_err,a_c_err]= DHI_profiles(r_ax(1:end-1),err_mean(1:end-1),Iz);

% figure(3)
% subplot(2,1,1)
% plot(r_ax,B_mean,'k','Linewidth',2)
% hold on
% errorbar(r_ax(1:2:end-1),B_prof(1:2:end-1),B_err(1:2:end-1),'ko','Linewidth',2)
% subplot(2,1,2)
% plot(r_ax,T_mean,'k','Linewidth',2)

% load dhi_map.mat
figure(4)
imagesc(x_twin,y_twin,den_num_full)
colormap jet; colorbar
set(gca,'Fontsize',14)
set(gca,'xtick',[0.043 0.048 0.052])


% mdsconnect('zappa.zap');
% mdsopen('zaphd',shotnum);
% m_0_p0 = mdsvalue('\m_0_p0');
% time = mdsvalue('dim_of(\m_0_p0)');
% time_holo = 45E-6;
% 
% 
% ind_time = find((time_holo - time)<0,1,'first');
% Bwall = m_0_p0(ind_time);
% rwall = 4*2.54/100;
% Iz = Bwall*2*pi*rwall/mu0;

[B_av,T_av,v_drift_av,den_char_av,B_char_av,...
    T_char_av,a_char_av] = DHI_profiles(r_ax,n_mean.*2E6,Iz);

[B_er,T_er,v_drift_er,den_char_er,B_char_er,...
    T_char_er,a_char_er] = DHI_profiles(r_ax,(n_mean(1:end-1) +err_mean).*2E6,Iz);

[B_d,T_d,v_drift_d,den_char_d,B_char_d,...
    T_char_d,a_char_d] = DHI_profiles(r_ax,(n_mean(1:end-1) -err_mean).*2E6,Iz);

erbar = abs(T_av-T_d);

figure(5)
plot(r_ax,T_av,'k','Linewidth',2)
hold on
errorbar(r_ax(1:5:end-1),T_av(1:5:end-1),erbar(1:5:end-1),'ko','LineWidth',2)

figure(6)
plot(r_ax,T_av)
hold on
plot(r_ax,T_er)
plot(r_ax,T_d)