clear all; clc;
close all;

directory_load = 'C:\Dropbox\MPROSS\DHIpaper\Coding_rebuttal\synthetic_holograms\profiles_lor_monteCarlo_nonoise.mat';
directory_save = 'C:\Dropbox\MPROSS\Dissertation\diss_latex_v2';

load(directory_load);
radius_NN{1} = radius;
den_int_NN{1} = den_int/10;
rad_rec_NN{1} = rad_rec;
den_int_rec_NN{1} = den_int_rec/10;
den_num_coarse_NN{1} = den_num_coarse/10;
den_num_rec_NN{1} = den_num_rec/10;
noise_image_NN{1} = noise_image;

%% computing uncertainty without noise:
% Lorentz:
error_N_NN{1} = den_int_rec_NN{1} - interp1(radius_NN{1},den_int_NN{1},rad_rec_NN{1});
error_N_NN{1}(find(isnan(error_N_NN{1}))) = 0;
error_n_NN{1} = abel_invert_uncertainty(rad_rec_NN{1},error_N_NN{1});

bar = 20;
fnt = 16;

fig1 = figure(1); hold on;
subplot(3,1,1); hold on;
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_tot_lor,'r');
% plot(linspace(0,0.011,length(den_num_coarse_NN{1})),den_num_coarse_NN{1},'b','linewidth',3);
errorbar(rad_rec_NN{1}(1:bar:end)*1000,den_num_rec_NN{1}(1:bar:end)/1e23,error_n_NN{1}(1:bar:end)/1e23,'r','linewidth',3);
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_mean_lor,'g');
xlim([0 11]);
ylim([0 2]);
% xlabel('radius [m]','fontsize',fnt);
ylabel('n_e [10^{23} m^{-3}]','fontsize',fnt);
% title('Lorentzian','fontsize',fnt);
set(gca,'fontsize',fnt);


a = 0.004;
n_0_1 = 0.75e23;
n_0_2 = 1e23;
n_0_3 = 1.25e23;


lorentz = (1/pi)*(a./(rad_rec_NN{1}.^2+a^2));
lorentz_zero = lorentz-min(lorentz);
lorentz_norm = lorentz_zero/max(lorentz_zero);
den_test{1} = n_0_1*lorentz_norm;
den_test{2} = n_0_2*lorentz_norm;
den_test{3} = n_0_3*lorentz_norm;
% den_test_2 = 
plt1 = plot(rad_rec_NN{1}*1000,den_test{1}/1e23,'color','g','linewidth',2);
plt2 = plot(rad_rec_NN{1}*1000,den_test{2}/1e23,'color','b','linestyle',':','linewidth',2);
plt3 = plot(rad_rec_NN{1}*1000,den_test{3}/1e23,'color','k','linestyle','--','linewidth',2);
leg1 = legend([plt1 plt2 plt3],'h = 1, n_0 = 0.75e23','h = 1, n_0 = 1e23','h = 1, n_0 = 1.25e23');
I = 1e5;

for i = 1:size(den_test,2)
    
    if size(den_test{i},1)>size(den_test{i},2)
    else
          den_test{i} = den_test{i}';
    end
    
    eval(['[B_prof_test{',num2str(i),'},T_prof_test{',num2str(i),...
        '},v_drift_test(',num2str(i),'),den_char_test(',num2str(i),...
        '),B_char_test(',num2str(i),'),T_char_test(',num2str(i),...
        '),a_char_test(',num2str(i),...
        ')] = DHI_profiles(rad_rec_NN{1},den_test{i},I);']);   
end

subplot(3,1,2); hold on;
plot(rad_rec_NN{1}*1000,B_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1}*1000,B_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1}*1000,B_prof_test{3},'color','k','linestyle','--','linewidth',2);
ylabel('B [T]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 11]);

subplot(3,1,3); hold on;
plot(rad_rec_NN{1}*1000,T_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{3},'color','k','linestyle','--','linewidth',2);
ylabel('T_i = T_e [eV]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlabel('radius [mm]','fontsize',fnt);
xlim([0 11]);

% saveas(fig1,[directory_save,'\DHI_equil_shape_1.png']);

return;

%%
clear all; clc;
% close all;

directory_load = 'C:\Dropbox\MPROSS\DHIpaper\Coding_rebuttal\synthetic_holograms\profiles_lor_monteCarlo_nonoise.mat';
directory_save = 'C:\Dropbox\MPROSS\Dissertation\diss_latex_v2';

load(directory_load);
radius_NN{1} = radius;
den_int_NN{1} = den_int/10;
rad_rec_NN{1} = rad_rec;
den_int_rec_NN{1} = den_int_rec/10;
den_num_coarse_NN{1} = den_num_coarse/10;
den_num_rec_NN{1} = den_num_rec/10;
noise_image_NN{1} = noise_image;

%% computing uncertainty without noise:
% Lorentz:
error_N_NN{1} = den_int_rec_NN{1} - interp1(radius_NN{1},den_int_NN{1},rad_rec_NN{1});
error_N_NN{1}(find(isnan(error_N_NN{1}))) = 0;
error_n_NN{1} = abel_invert_uncertainty(rad_rec_NN{1},error_N_NN{1});

a = 0.004;
n_0 = 1e23;
h = [1/2 1 2];

for j = 1:length(h)
    lorentz = (1/pi)*(a./(h(j)*rad_rec_NN{1}.^2+a^2));
    lorentz_zero = lorentz-min(lorentz);
    lorentz_norm = lorentz_zero/max(lorentz_zero);
    den_test{j} = n_0*lorentz_norm;    
end

bar = 20;
fnt = 16;

fig2 = figure(2); hold on;
subplot(3,1,1); hold on;
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_tot_lor,'r');
% plot(linspace(0,0.011,length(den_num_coarse_NN{1})),den_num_coarse_NN{1},'b','linewidth',3);
errorbar(rad_rec_NN{1}(1:bar:end)*1000,den_num_rec_NN{1}(1:bar:end)/1e23,error_n_NN{1}(1:bar:end)/1e23,'r','linewidth',3);
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_mean_lor,'g');
xlim([0 11]);
ylim([0 2]);
% xlabel('radius [m]','fontsize',fnt);
ylabel('n_e [10^{23} m^{-3}]','fontsize',fnt);
title('Lorentzian','fontsize',fnt);
set(gca,'fontsize',fnt);

plt1 = plot(rad_rec_NN{1}*1000,den_test{1}/1e23,'color','g','linewidth',2);
plt2 = plot(rad_rec_NN{1}*1000,den_test{2}/1e23,'color','b','linestyle',':','linewidth',2);
plt3 = plot(rad_rec_NN{1}*1000,den_test{3}/1e23,'color','k','linestyle','--','linewidth',2);
leg1 = legend([plt1 plt2 plt3],'h = 0.5, n_0 = 1e23','h = 1, n_0 = 1e23','h = 2, n_0 = 1e23');

I = 1e5;

for i = 1:size(den_test,2)
    
    if size(den_test{i},1)>size(den_test{i},2)
    else
          den_test{i} = den_test{i}';
    end
    eval(['[B_prof_test{',num2str(i),'},T_prof_test{',num2str(i),...
        '},v_drift_test(',num2str(i),'),den_char_test(',num2str(i),...
        '),B_char_test(',num2str(i),'),T_char_test(',num2str(i),...
        '),a_char_test(',num2str(i),...
        ')] = DHI_profiles(rad_rec_NN{1},den_test{i},I);']);   
end

subplot(3,1,2); hold on;
plot(rad_rec_NN{1}*1000,B_prof_test{1},'color','g','linewidth',2);
% scatter((a_char_test(1)),B_char_test(1));
plot(rad_rec_NN{1}*1000,B_prof_test{2},'color','b','linestyle',':','linewidth',2);
% scatter((a_char_test(2)),B_char_test(2));
plot(rad_rec_NN{1}*1000,B_prof_test{3},'color','k','linestyle','--','linewidth',2);
% scatter((a_char_test(3)),B_char_test(3));
ylabel('B [T]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 11]);

subplot(3,1,3); hold on;
plot(rad_rec_NN{1}*1000,T_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{3},'color','k','linestyle','--','linewidth',2);
ylabel('T_i = T_e [eV]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 11]);
xlabel('radius [mm]','fontsize',fnt);

% saveas(fig2,[directory_save,'\DHI_equil_shape_2.png']);

%%
clear all; clc;
% close all;

directory_load = 'C:\Dropbox\MPROSS\DHIpaper\Coding_rebuttal\synthetic_holograms\profiles_lor_monteCarlo_nonoise.mat';
directory_save = 'C:\Dropbox\MPROSS\Dissertation\diss_latex_v2';

load(directory_load);
radius_NN{1} = radius;
den_int_NN{1} = den_int;
rad_rec_NN{1} = rad_rec;
den_int_rec_NN{1} = den_int_rec;
den_num_coarse_NN{1} = den_num_coarse;
den_num_rec_NN{1} = den_num_rec;
noise_image_NN{1} = noise_image;

%% computing uncertainty without noise:
% Lorentz:
error_N_NN{1} = den_int_rec_NN{1} - interp1(radius_NN{1},den_int_NN{1},rad_rec_NN{1});
error_N_NN{1}(find(isnan(error_N_NN{1}))) = 0;
error_n_NN{1} = abel_invert_uncertainty(rad_rec_NN{1},error_N_NN{1});

a = 0.004;
n_0 = linspace(.75e24,1.25e24,25);%[.5 .75 1 1.25 1.5]*1e24;
h = linspace(0.5,2,25);%[1/2 1 2];

for i = 1:length(n_0)
    for j = 1:length(h)
        lorentz = (1/pi)*(a./(h(j)*rad_rec_NN{1}.^2+a^2));
        lorentz_zero = lorentz-min(lorentz);
        lorentz_norm = lorentz_zero/max(lorentz_zero);
        den_test{i,j} = n_0(i)*lorentz_norm;
        maxes(i,j) = max(den_test{i,j});
    end
end

bar = 20;
fnt = 16;

% fig3 = figure(3); hold on;
% subplot(1,3,1); hold on;
% pcolor(n_0,h,den_test{i,j});shading interp;
% 
% % subplot(3,1,1); hold on;
% % errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_tot_lor,'r');
% % plot(linspace(0,0.011,length(den_num_coarse_NN{1})),den_num_coarse_NN{1},'b','linewidth',3);
% errorbar(rad_rec_NN{1}(1:bar:end),den_num_rec_NN{1}(1:bar:end),error_n_NN{1}(1:bar:end),'r','linewidth',3);
% % errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_mean_lor,'g');
% xlim([0 0.011]);
% ylim([0 2e24]);
% % xlabel('radius [m]','fontsize',fnt);
% ylabel('n_e [m^{-3}]','fontsize',fnt);
% title('Lorentzian','fontsize',fnt);
% set(gca,'fontsize',fnt);
% xlabel('radius [m]','fontsize',fnt);
% 
% plot(rad_rec_NN{1},den_test{1},'color','g','linewidth',2);
% plot(rad_rec_NN{1},den_test{2},'color','b','linestyle',':','linewidth',2);
% plot(rad_rec_NN{1},den_test{3},'color','k','linestyle','--','linewidth',2);

I = 1e5;

for i = 1:length(n_0)
    for j = 1:length(h)
        eval(['[B_prof_test{',num2str(i),',',num2str(j),'},T_prof_test{',num2str(i),',',num2str(j),...
            '},v_drift_test(',num2str(i),',',num2str(j),'),den_char_test(',num2str(i),',',num2str(j),...
            '),B_char_test(',num2str(i),',',num2str(j),'),T_char_test(',num2str(i),',',num2str(j),...
            '),a_char_test(',num2str(i),',',num2str(j),...
            ')] = DHI_profiles(rad_rec_NN{1},den_test{i,j},I);']);
    end
end

figure(11); hold on;
pcolor(h,n_0,den_char_test); shading interp;
colorbar;
colormap jet;

figure(12); hold on;
pcolor(h,n_0,B_char_test); shading interp;
colorbar;
colormap jet;

figure(13); hold on;
pcolor(h,n_0,T_char_test); shading interp;
colorbar;
colormap jet;

subplot(3,1,2); hold on;
plot(rad_rec_NN{1},B_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1},B_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1},B_prof_test{3},'color','k','linestyle','--','linewidth',2);
ylabel('B [T]','fontsize',fnt);

subplot(3,1,3); hold on;
plot(rad_rec_NN{1},T_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1},T_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1},T_prof_test{3},'color','k','linestyle','--','linewidth',2);
ylabel('T_i = T_e [eV]','fontsize',fnt);

% saveas(fig3,[directory_save,'\DHI_equil_shape_3.png']);




%%
clear all; clc;
% close all;

directory_load = 'C:\Dropbox\MPROSS\DHIpaper\Coding_rebuttal\synthetic_holograms\profiles_hyp_monteCarlo_nonoise.mat';
directory_save = 'C:\Dropbox\MPROSS\Dissertation\diss_latex_v2';

load(directory_load);
radius_NN{1} = radius;
den_int_NN{1} = den_int/10;
rad_rec_NN{1} = rad_rec;
den_int_rec_NN{1} = den_int_rec/10;
den_num_coarse_NN{1} = den_num_coarse/10;
den_num_rec_NN{1} = den_num_rec/10;
noise_image_NN{1} = noise_image;

%% computing uncertainty without noise:
% Lorentz:
error_N_NN{1} = den_int_rec_NN{1} - interp1(radius_NN{1},den_int_NN{1},rad_rec_NN{1});
error_N_NN{1}(find(isnan(error_N_NN{1}))) = 0;
error_n_NN{1} = abel_invert_uncertainty(rad_rec_NN{1},error_N_NN{1});

a = 0.0055;
n_0 = 1e23;
h = [5 10 75 10];

for j = 1:length(h)
        
%     magn = 10;
%     tanh_prof = (1-tanh(h(j)*(rad_rec_NN{1}-a)/a));
%     tanh_prof_zero = tanh_prof-min(tanh_prof);
%     tanh_prof_norm = tanh_prof_zero/max(tanh_prof_zero);
%     den_test{j} = n_0*(tanh_prof_zero/max(tanh_prof_zero));
    if j <= 3
        den_test{j} = n_0*0.5*(1-tanh(h(j)*(rad_rec_NN{1}-a)/a));
    else
        a = .0045;
        den_test{j} = n_0*0.5*(1-tanh(h(j)*(rad_rec_NN{1}-a)/a));
    end
end

bar = 20;
fnt = 16;

fig3 = figure(3); hold on;
subplot(3,1,1); hold on;
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_tot_lor,'r');
% plot(linspace(0,0.011,length(den_num_coarse_NN{1})),den_num_coarse_NN{1},'b','linewidth',3);
errorbar(rad_rec_NN{1}(1:bar:end)*1000,den_num_rec_NN{1}(1:bar:end)/1e23,error_n_NN{1}(1:bar:end)/1e23,'r','linewidth',3);
% errorbar(rad_rec_NN{1},den_num_rec_NN{1},error_mean_lor,'g');
xlim([0 11]);
ylim([0 2]);
% xlabel('radius [m]','fontsize',fnt);
ylabel('n_e [10^{23} m^{-3}]','fontsize',fnt);
% title('Hyperbolic Tangent','fontsize',fnt);
set(gca,'fontsize',fnt);

plt1 = plot(rad_rec_NN{1}*1000,den_test{1}/1e23,'color','g','linewidth',2);
plt2 = plot(rad_rec_NN{1}*1000,den_test{2}/1e23,'color','b','linestyle',':','linewidth',2);
plt3 = plot(rad_rec_NN{1}*1000,den_test{3}/1e23,'color','k','linestyle','--','linewidth',2);
plt4 = plot(rad_rec_NN{1}*1000,den_test{4}/1e23,'color','c','linestyle','--','linewidth',2);
leg1 = legend([plt1 plt2 plt3 plt4],'k = 5, a = 0.0055','h = 10, a = 0.0055','h = 75, a = 0.0055','h = 10, a = 0.0045');

I = 1e5;

for i = 1:size(den_test,2)
    
    den_test{i} = (den_test{i}');
    
    eval(['[B_prof_test{',num2str(i),'},T_prof_test{',num2str(i),...
        '},v_drift_test(',num2str(i),'),den_char_test(',num2str(i),...
        '),B_char_test(',num2str(i),'),T_char_test(',num2str(i),...
        '),a_char_test(',num2str(i),...
        ')] = DHI_profiles(rad_rec_NN{1},den_test{i},I);']);   
end

subplot(3,1,2); hold on;
plot(rad_rec_NN{1}*1000,B_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1}*1000,B_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1}*1000,B_prof_test{3},'color','k','linestyle','--','linewidth',2);
plot(rad_rec_NN{1}*1000,B_prof_test{4},'color','c','linestyle','--','linewidth',2);
ylabel('B [T]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 11]);



subplot(3,1,3); hold on;
plot(rad_rec_NN{1}*1000,T_prof_test{1},'color','g','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{2},'color','b','linestyle',':','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{3},'color','k','linestyle','--','linewidth',2);
plot(rad_rec_NN{1}*1000,T_prof_test{4},'color','c','linestyle','--','linewidth',2);
ylabel('T_i = T_e [eV]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 11]);
xlabel('radius [mm]','fontsize',fnt);


% saveas(fig3,[directory_save,'\DHI_equil_shape_3.png']);