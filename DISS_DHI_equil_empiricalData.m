clear all; clc;
close all;

shotnum = 170504050;

if shotnum == 160524021
    load(['M:\Users\Michael_Ross\Coding\Holography_Reconstruction\dhi_error_data.mat']);
else
    load(['M:\Users\Michael_Ross\Coding\Holography_Reconstruction\dhi_error_data_',num2str(shotnum),'.mat']);
end

r_w = 4*2.54/100;
mu_0 = (4*pi)*1e-7;
k_boltz =1.38e-23; % J/K
m_i = 1.67e-27;

[time,m_0_p0] = MDS_MATLAB_python(shotnum,'m_0_p0');
[time_holo] = MDS_MATLAB_python_digitizers(shotnum,...
    'ZAPHD::TOP.DIGITIZERSHD.L4222:ZAP_4222_B:CHANNEL_1');
ind_time = find((time_holo - time)<0,1,'first');
B_w = m_0_p0(ind_time);
I = B_w*2*pi*r_w/mu_0;

k =2;
right = 1;
if right == 1
rad_r{k} = rad_l{k};
den_num_axial_r_plt{k} = den_num_axial_l_plt{k};
n_error_out_interp_r{k} = n_error_out_interp_l{k};
end

bar_spacing = 7;

[B_prof,T_prof,v_drift,den_char,B_char,...
    T_char,a_char] = DHI_profiles(rad_r{k},den_num_axial_r_plt{k},I);

fnt = 16;
figure(1); hold on;
subplot(3,1,1); hold on;
plot(rad_r{k}*1000,den_num_axial_r_plt{k}/1e23,'color','k','linewidth',3);
errorbar(rad_r{k}(1:bar_spacing:end)*1000,den_num_axial_r_plt{k}(1:bar_spacing:end)/1e23,n_error_out_interp_r{k}(1:bar_spacing:end)/1e23,'color','k','linestyle','none','linewidth',3);
ylabel('n_e [10^{23} m^{-3}]','fontsize',fnt);
sct1 = scatter(0,den_char/1e23,'r');
line([a_char a_char]*1000,[0 4],'color','r','linestyle',':','linewidth',2);
xlim([0 11]);
ylim([0 4]);
% title(['right = ',num2str(right),', k = ',num2str(k),', shotnum = ',num2str(shotnum)]);
tit1 = title(num2str(shotnum),'fontsize',8);
set(gca,'fontsize',fnt);
set(tit1,'fontsize',8);

subplot(3,1,2); hold on;
plot(rad_r{k}*1000,B_prof,'color','k','linewidth',3);
ylabel('B [T]','fontsize',fnt);
sct1 = scatter(a_char*1000,B_char,'r');
line([a_char a_char]*1000,[0 ceil(B_char)],'color','r','linestyle',':','linewidth',2);
xlim([0 11]);
ylim([0 ceil(B_char)]);
set(gca,'fontsize',fnt);

subplot(3,1,3); hold on;
plot(rad_r{k}*1000,T_prof,'color','k','linewidth',3);
ylabel('T_e [eV]','fontsize',fnt);
sct1 = scatter(0,T_char,'r');
line([a_char a_char]*1000,[0 ceil(T_char)],'color','r','linestyle',':','linewidth',2);
xlim([0 11]);
ylim([0 ceil(T_char)]);
set(gca,'fontsize',fnt);
xlabel('radius [mm]','fontsize',fnt);