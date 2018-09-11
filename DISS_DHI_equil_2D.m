% clear all; clc;
% close all;

shotnum = 171005023;%160524018;

% dhi_inverted = MDS_MATLAB_python_dhihd(shotnum,'dhi_inverted');
% centroid_abs = MDS_MATLAB_python_dhihd(shotnum,'centroid_abs');
% x_twin = MDS_MATLAB_python_dhihd(shotnum,'x_twin');
% y_twin = MDS_MATLAB_python_dhihd(shotnum,'y_twin');
x_twin = x_twin_red_adj;
y_twin = y_twin_red_adj;
dhi_inverted = den_num_full;

mu0 = (4*pi)*1e-7;
% [time,m_0_p0] = MDS_MATLAB_python(shotnum,'m_0_p0');
% [time_holo] = MDS_MATLAB_python_digitizers(shotnum,...
%     'ZAPHD::TOP.DIGITIZERSHD.L4222:ZAP_4222_B:CHANNEL_1');
mdsconnect('zappa.zap');
mdsopen('zaphd',shotnum);
m_0_p0 = mdsvalue('\m_0_p0');
time = mdsvalue('dim_of(\m_0_p0)');
time_holo = 45E-6;


ind_time = find((time_holo - time)<0,1,'first');
Bwall = m_0_p0(ind_time);
rwall = 4*2.54/100;
Iz = Bwall*2*pi*rwall/mu0;

del_y = y_twin(2) - y_twin(1);

for i = 1:length(x_twin)
    r_l{i} = 0:del_y:(centroid_abs(i)-1)*del_y;
    n_l{i} = flipud(dhi_inverted(1:centroid_abs(i),i));
    r_r{i} = 0:del_y:(length(y_twin)-centroid_abs(i)-1)*del_y;
    n_r{i} = dhi_inverted(centroid_abs(i)+1:length(y_twin),i);
end

for i = 1:length(x_twin)
    
    if size(n_l{i},1)>size(n_l{i},2)
    else
          n_l{i} = n_l{i}';
    end
    if size(n_r{i},1)>size(n_r{i},2)
    else
          n_r{i} = n_r{i}';
    end
    [B_prof_l{i},T_prof_l{i},v_drift_l(i),den_char_l(i),B_char_l(i),...
    T_char_l(i),a_char_l(i)] = DHI_profiles(r_l{i},n_l{i},Iz);

    [B_prof_r{i},T_prof_r{i},v_drift_r(i),den_char_r(i),B_char_r(i),...
    T_char_r(i),a_char_r(i)] = DHI_profiles(r_r{i},n_r{i},Iz);
end

B_mat = zeros(size(dhi_inverted));
T_mat = zeros(size(dhi_inverted));
for i = 1:length(x_twin)
    B_mat(1:length(B_prof_l{i}),i) = flipud(B_prof_l{i});
%     a_l(i) = find(ca_char_l(i)-r_l{i}<0,1,'first');
    B_mat(length(B_prof_l{i})+1:end,i) = (B_prof_r{i});
%     a_r(i) = find(a_char_r(i)-r_r{i}<0,1,'first');
    T_mat(1:length(T_prof_l{i}),i) = flipud(T_prof_l{i});
    T_mat(length(T_prof_l{i})+1:end,i) = (T_prof_r{i});
end

for i = 1:length(x_twin)
    a_l(i) = y_twin(centroid_abs(i))-a_char_l(i);
    a_r(i) = y_twin(centroid_abs(i))+a_char_r(i);
end

fnt = 24;

ccol = 30;
%%
fig1 = figure(1); hold on;
imagesc(x_twin,y_twin,dhi_inverted);
h32 = scatter(x_twin,y_twin(centroid_abs),'k.');
set(h32,'sizedata',1000);
sct1 = scatter(x_twin,a_l,'m.');
set(sct1,'sizedata',1000);
sct2 = scatter(x_twin,a_r,'m.');
set(sct2,'sizedata',1000);
set(gca,'ydir','normal');
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
title('Electron number density, n_e [m^{-3}]','fontsize',fnt);
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .080 .083]);

figure(2); hold on;
plot(r_l{ccol},n_l{ccol});
plot(r_r{ccol},n_r{ccol},'r');

figure(3); hold on;
waterfall(B_mat');
colormap jet;

figure(4); hold on;
imagesc(x_twin,y_twin,B_mat);
% h32 = scatter(x_twin,y_twin(centroid_abs),'k','.');
% set(h32,'sizedata',100);
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
title('Magnetic Field [T]','fontsize',fnt);
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .080 .083]);
% set(gca,'clim',[0 15]);

figure(5); hold on;
waterfall(T_mat');
colormap jet;

figure(6); hold on;
imagesc(x_twin,y_twin,T_mat);
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
title('T_e = T_i [eV]','fontsize',fnt);
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .080 .083]);
% set(gca,'clim',[0 5000]);
% save(['M:\Users\Michael_Ross\Coding\video_read\computed',num2str(shotnum),'.mat'],'B_mat','T_mat');


fnt = 16;
figure(7);
subplot(3,1,1); hold on;
plot(r_r{ccol},n_r{ccol},'color','k','linewidth',3);
ylabel('n_e [m^{-3}]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 0.012]);

subplot(3,1,2); hold on;
plot(r_r{ccol},B_prof_r{ccol},'color','k','linewidth',3);
ylabel('B [T]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 0.012]);

subplot(3,1,3); hold on;
plot(r_r{ccol},T_prof_r{ccol},'color','k','linewidth',3);
ylabel('T_i = T_e [eV]','fontsize',fnt);
set(gca,'fontsize',fnt);
xlim([0 0.012]);


a2_vec = linspace(0,0.015,100);
fnt = 18;
scatter_colors = colormap(colorcube);
%% density vs. size:
fig10 = figure(10); hold on;

n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2)*(1);
% n2_a2_vec_N2 = n1*(a2_vec/a1).^(-2)*(1/4);
% n2_a2_vec_N3 = n1*(a2_vec/a1).^(-2)*(1/8);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter(a_char_l(i),den_char_l(i),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),den_char_r(i),'k','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end

p1 = patch([0.005 0.005 0.015 0.015],[1e22 1e23 1e23 1e22],'r','EdgeColor','r');
alpha(p1,0.25);

% plt1 = plot(a2_vec,n2_a2_vec,'k');
plt2 = plot(a2_vec,n2_a2_vec_N1,'r');
% plt3 = plot(a2_vec,n2_a2_vec_N2,'b');
% plt4 = plot(a2_vec,n2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak density [m^{-3}]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 6e23]);

% leg1 = legend([plt1 plt2 plt3 plt4],'N_2/N_1 = 1','N_2/N_1 = 1/2','N_2/N_1 = 1/4','N_2/N_1 = 1/8');


gamma = 5/3;
%% field vs. size:
fig11 = figure(11); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    sct1 = scatter(a_char_l(i),B_char_l(i),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),B_char_r(i),'k','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end

p1 = patch([0.005 0.005 0.015 0.015],[.3 1.8 1.8 0.3],'r','EdgeColor','r');
alpha(p1,0.25);
% plot(a2_vec,B2_a2_vec,'k');
plot(a2_vec,B2_a2_vec_N1,'r');
% plot(a2_vec,B2_a2_vec_N2,'b');
% plot(a2_vec,B2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak field [T]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 20]);

% saveas(fig11,[directory_load,'\DHI_confinement_11.png']);



%% temperature vs. size:
fig12 = figure(12); hold on;

T2_a2_vec_N1 = T_char_l(i)*(a2_vec/a_char_l(i)).^(-2*(gamma-1))*(1)^(gamma-1);
% T2_a2_vec_N2 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/4)^(gamma-1);
% T2_a2_vec_N3 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/8)^(gamma-1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter(a_char_l(i),T_char_l(i),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),T_char_r(i),'r','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end
p2 = patch([0.005 0.005 0.015 0.015],[20 220 220 20],'r','EdgeColor','r');
alpha(p2,0.25);
% plot(a2_vec,T2_a2_vec,'k');
plot(a2_vec,T2_a2_vec_N1,'r');
% plot(a2_vec,T2_a2_vec_N2,'b');
% plot(a2_vec,T2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak temperature [eV]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 7500]);

% saveas(fig12,[directory_load,'\DHI_confinement_12.png']);


%% Computing linear density:
for i = length(x_twin):-1:1
    
    % computing limits of integration (only want to integrate density out
    % to a:
    int_upper_l(i) = find((r_l{i}-a_char_l(i))>0,1,'first');
    int_upper_r(i) = find((r_r{i}-a_char_r(i))>0,1,'first');

    del_y = y_twin(2)-y_twin(1);
    
    rad_linear_l{i} = 0:del_y:(length(n_l{i})-1)*del_y;  % del_y in empirical reconstruction coords.
    rad_linear_r{i} = 0:del_y:(length(n_r{i})-1)*del_y; % del_y in empirical reconstruction coords.
    
    N_linear_density_l_tmp = cumtrapz(rad_linear_l{i},2*pi*n_l{i}.*rad_linear_l{i}');
    N_linear_density_l(i) = N_linear_density_l_tmp(int_upper_l(i));
    
    N_linear_density_r_tmp = cumtrapz(rad_linear_r{i},2*pi*n_r{i}.*rad_linear_r{i}');
    N_linear_density_r(i) = N_linear_density_r_tmp(int_upper_r(i));
end


fig19 = figure(19); hold on;

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter(a_char_l(i),N_linear_density_l(i),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),N_linear_density_r(i),'k','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end

xlabel('characteristic size [A]','fontsize',fnt);
ylabel('linear density [m^{-1}]','fontsize',fnt);
set(gca,'fontsize',fnt);

% saveas(fig19,[directory_load,'\DHI_confinement_19.png']);




%% density vs. size adjusted for linear density:
fig20 = figure(20); hold on;

n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2);
% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2)*(1);
% n2_a2_vec_N2 = n1*(a2_vec/a1).^(-2)*(1/4);
% n2_a2_vec_N3 = n1*(a2_vec/a1).^(-2)*(1/8);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter(a_char_l(i),den_char_l(i)*N_linear_density_l(1)/N_linear_density_l(i),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),den_char_r(i)*N_linear_density_l(1)/N_linear_density_r(i),'k','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end

p1 = patch([0.005 0.005 0.015 0.015],[1e22 1e23 1e23 1e22],'r','EdgeColor','r');
alpha(p1,0.25);

% plt1 = plot(a2_vec,n2_a2_vec,'k');
plt2 = plot(a2_vec,n2_a2_vec_N1,'r');
% plt3 = plot(a2_vec,n2_a2_vec_N2,'b');
% plt4 = plot(a2_vec,n2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak density [m^{-3}]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 6e23]);

% leg1 = legend([plt1 plt2 plt3 plt4],'N_2/N_1 = 1','N_2/N_1 = 1/2','N_2/N_1 = 1/4','N_2/N_1 = 1/8');


gamma = 5/3;
%% field vs. size:
fig21 = figure(21); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    sct1 = scatter(a_char_l(i),B_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(exponent_1),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
    sct2 = scatter(a_char_r(i),B_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(exponent_1),'k','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end

p1 = patch([0.005 0.005 0.015 0.015],[.3 1.8 1.8 0.3],'r','EdgeColor','r');
alpha(p1,0.25);
% plot(a2_vec,B2_a2_vec,'k');
plot(a2_vec,B2_a2_vec_N1,'r');
% plot(a2_vec,B2_a2_vec_N2,'b');
% plot(a2_vec,B2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak field [T]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 20]);

% saveas(fig21,[directory_load,'\DHI_confinement_21.png']);


%% temperature vs. size:
fig22 = figure(22); hold on;

T2_a2_vec_N1 = T_char_l(i)*(a2_vec/a_char_l(i)).^(-2*(gamma-1))*(1)^(gamma-1);
% T2_a2_vec_N2 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/4)^(gamma-1);
% T2_a2_vec_N3 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/8)^(gamma-1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
%     sct1 = scatter(a_char_l(i),T_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(gamma-1),'k','d','filled');
    sct1 = scatter(a_char_l(i),T_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(gamma-1),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     sct2 = scatter(a_char_r(i),T_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(gamma-1),'r','s','filled');
    sct2 = scatter(a_char_r(i),T_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(gamma-1),'r','s','filled');
    if i < size(scatter_colors,1)
        set(sct2,'CData',scatter_colors(i,:));
    else
        set(sct2,'CData',scatter_colors(end,:));
    end
    set(sct1,'MarkerEdgeColor','k');
    set(sct2,'MarkerEdgeColor','k');
end
p2 = patch([0.005 0.005 0.015 0.015],[20 220 220 20],'r','EdgeColor','r');
alpha(p2,0.25);
% plot(a2_vec,T2_a2_vec,'k');
plot(a2_vec,T2_a2_vec_N1,'r');
% plot(a2_vec,T2_a2_vec_N2,'b');
% plot(a2_vec,T2_a2_vec_N3,'m');

xlabel('characteristic size [m]','fontsize',fnt);
ylabel('peak temperature [eV]','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([0 0.015]);
ylim([0 7500]);

% saveas(fig12,[directory_load,'\DHI_confinement_12.png']);



mar_size = 200;
x_line = linspace(0,4,100);
y_line = linspace(0,4,100);

%% density vs. size adjusted for linear density:
fig23 = figure(23); hold on;

% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(1)).^(-2);
% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2)*(1);
% n2_a2_vec_N2 = n1*(a2_vec/a1).^(-2)*(1/4);
% n2_a2_vec_N3 = n1*(a2_vec/a1).^(-2)*(1/8);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2),(den_char_l(i)/den_char_l(1))...
        *(N_linear_density_l(i)/N_linear_density_l(1))^(-1),'k','d','filled');
%     if i < size(scatter_colors,1)
%         set(sct1,'CData',scatter_colors(i,:));
%     else
%         set(sct1,'CData',scatter_colors(end,:));
%     end
    set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2),(den_char_r(i)/den_char_l(1))...
%         *(N_linear_density_r(i)/N_linear_density_l(1))^(-1),'k','s','filled');    
%     if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

plot(x_line,y_line,'k:');
% p1 = patch([0.005 0.005 0.015 0.015],[1e22 1e23 1e23 1e22],'r','EdgeColor','r');
% alpha(p1,0.25);

% plt1 = plot(a2_vec,n2_a2_vec,'k');
% plt2 = plot(a2_vec,n2_a2_vec_N1,'r');
% plt3 = plot(a2_vec,n2_a2_vec_N2,'b');
% plt4 = plot(a2_vec,n2_a2_vec_N3,'m');

xlabel('(a_2/a_1)^{-2}','fontsize',fnt);
ylabel('(n_2/n_1) (N_2/N_1) ^{-1}','fontsize',fnt);
title('density scaling','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
lim_max = 2;
lim_min = 0.5;
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% leg1 = legend([plt1 plt2 plt3 plt4],'N_2/N_1 = 1','N_2/N_1 = 1/2','N_2/N_1 = 1/4','N_2/N_1 = 1/8');


gamma = 5/3;
%% field vs. size:
fig24 = figure(24); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
% B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-gamma),(B_char_l(i)/B_char_l(1))*...
        (N_linear_density_l(i)/N_linear_density_l(1))^(-exponent_1),'k','d','filled');
%     if i < size(scatter_colors,1)
%         set(sct1,'CData',scatter_colors(i,:));
%     else
%         set(sct1,'CData',scatter_colors(end,:));
%     end
    set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-gamma),(B_char_r(i)/B_char_l(1))*...
%         (N_linear_density_r(i)/N_linear_density_l(1))^(-exponent_1),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

% p1 = patch([0.005 0.005 0.015 0.015],[.3 1.8 1.8 0.3],'r','EdgeColor','r');
% alpha(p1,0.25);
% % plot(a2_vec,B2_a2_vec,'k');
% plot(a2_vec,B2_a2_vec_N1,'r');
% plot(a2_vec,B2_a2_vec_N2,'b');
% plot(a2_vec,B2_a2_vec_N3,'m');
plot(x_line,y_line,'k:');


xlabel('(a_2/a_1)^{-\gamma}','fontsize',fnt);
ylabel('(B_2/B_1) (N_2/N_1) ^{-(\gamma^2-\gamma)/(2(\gamma-1))}','fontsize',fnt);
title('Magnetic field scaling','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
lim_max = 2;
lim_min = 0.5;
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig21,[directory_load,'\DHI_confinement_21.png']);


%% temperature vs. size:
fig25 = figure(25); hold on;

% T2_a2_vec_N1 = T_char_l(i)*(a2_vec/a_char_l(i)).^(-2*(gamma-1))*(1)^(gamma-1);
% T2_a2_vec_N2 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/4)^(gamma-1);
% T2_a2_vec_N3 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/8)^(gamma-1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
%     sct1 = scatter(a_char_l(i),T_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(gamma-1),'k','d','filled');
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2*(gamma-1)),...
        (T_char_l(i)/T_char_l(1))*(N_linear_density_l(i)/N_linear_density_l(1))^(1-gamma),'k','d','filled');
%     if i < size(scatter_colors,1)
%         set(sct1,'CData',scatter_colors(i,:));
%     else
%         set(sct1,'CData',scatter_colors(end,:));
%     end
    set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter(a_char_r(i),T_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(gamma-1),'r','s','filled');
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2*(gamma-1)),...
%         (T_char_r(i)/T_char_l(1))*(N_linear_density_r(i)/N_linear_density_l(1))^(1-gamma),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end
% p2 = patch([0.005 0.005 0.015 0.015],[20 220 220 20],'r','EdgeColor','r');
% alpha(p2,0.25);
% plot(a2_vec,T2_a2_vec,'k');
% plot(a2_vec,T2_a2_vec_N1,'r');
% plot(a2_vec,T2_a2_vec_N2,'b');
% plot(a2_vec,T2_a2_vec_N3,'m');
plot(x_line,y_line,'k:');

xlabel('(a_2/a_1)^{-2(\gamma-1)}','fontsize',fnt);
ylabel('(T_2/T_1) (N_2/N_1)^{1-\gamma}','fontsize',fnt);
title('Temperature scaling','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);

lim_max = 2;
lim_min = 0.5;
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig12,[directory_load,'\DHI_confinement_12.png']);


%% subplots of normalized scaling:
scatter_colors = colormap(copper);
fig30 = figure(30);

mar_size = 200;
x_line = linspace(0,4,100);
y_line = linspace(0,4,100);

%% density vs. size adjusted for linear density:
subplot(1,3,1); hold on;

% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(1)).^(-2);
% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2)*(1);
% n2_a2_vec_N2 = n1*(a2_vec/a1).^(-2)*(1/4);
% n2_a2_vec_N3 = n1*(a2_vec/a1).^(-2)*(1/8);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2),(den_char_l(i)/den_char_l(1))...
        *(N_linear_density_l(i)/N_linear_density_l(1))^(-1),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2),(den_char_r(i)/den_char_l(1))...
%         *(N_linear_density_r(i)/N_linear_density_l(1))^(-1),'k','s','filled');    
%     if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

plot(x_line,y_line,'k:','linewidth',2);
% p1 = patch([0.005 0.005 0.015 0.015],[1e22 1e23 1e23 1e22],'r','EdgeColor','r');
% alpha(p1,0.25);

% plt1 = plot(a2_vec,n2_a2_vec,'k');
% plt2 = plot(a2_vec,n2_a2_vec_N1,'r');
% plt3 = plot(a2_vec,n2_a2_vec_N2,'b');
% plt4 = plot(a2_vec,n2_a2_vec_N3,'m');

xlabel('(a_2/a_1)^{-2}','fontsize',fnt);
ylabel('(n_2/n_1) (N_2/N_1) ^{-1}','fontsize',fnt);
title('Density','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
lim_max = 3.5;
lim_min = .5;
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% leg1 = legend([plt1 plt2 plt3 plt4],'N_2/N_1 = 1','N_2/N_1 = 1/2','N_2/N_1 = 1/4','N_2/N_1 = 1/8');


gamma = 5/3;
%% field vs. size:
subplot(1,3,2); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
% B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-gamma),(B_char_l(i)/B_char_l(1))*...
        (N_linear_density_l(i)/N_linear_density_l(1))^(-exponent_1),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-gamma),(B_char_r(i)/B_char_l(1))*...
%         (N_linear_density_r(i)/N_linear_density_l(1))^(-exponent_1),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

% p1 = patch([0.005 0.005 0.015 0.015],[.3 1.8 1.8 0.3],'r','EdgeColor','r');
% alpha(p1,0.25);
% % plot(a2_vec,B2_a2_vec,'k');
% plot(a2_vec,B2_a2_vec_N1,'r');
% plot(a2_vec,B2_a2_vec_N2,'b');
% plot(a2_vec,B2_a2_vec_N3,'m');
plot(x_line,y_line,'k:','linewidth',2);


xlabel('(a_2/a_1)^{-\gamma}','fontsize',fnt);
ylabel('(B_2/B_1) (N_2/N_1) ^{-\gamma^2/2}','fontsize',fnt);
title('Magnetic field','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig21,[directory_load,'\DHI_confinement_21.png']);


%% temperature vs. size:
subplot(1,3,3); hold on;

% T2_a2_vec_N1 = T_char_l(i)*(a2_vec/a_char_l(i)).^(-2*(gamma-1))*(1)^(gamma-1);
% T2_a2_vec_N2 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/4)^(gamma-1);
% T2_a2_vec_N3 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/8)^(gamma-1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
%     sct1 = scatter(a_char_l(i),T_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(gamma-1),'k','d','filled');
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2*(gamma-1)),...
        (T_char_l(i)/T_char_l(1))*(N_linear_density_l(i)/N_linear_density_l(1))^(1-gamma),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter(a_char_r(i),T_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(gamma-1),'r','s','filled');
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2*(gamma-1)),...
%         (T_char_r(i)/T_char_l(1))*(N_linear_density_r(i)/N_linear_density_l(1))^(1-gamma),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end
% p2 = patch([0.005 0.005 0.015 0.015],[20 220 220 20],'r','EdgeColor','r');
% alpha(p2,0.25);
% plot(a2_vec,T2_a2_vec,'k');
% plot(a2_vec,T2_a2_vec_N1,'r');
% plot(a2_vec,T2_a2_vec_N2,'b');
% plot(a2_vec,T2_a2_vec_N3,'m');
plot(x_line,y_line,'k:','linewidth',2);

xlabel('(a_2/a_1)^{-2(\gamma-1)}','fontsize',fnt);
ylabel('(T_2/T_1) (N_2/N_1)^{1-\gamma}','fontsize',fnt);
title('Temperature','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);


xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig12,[directory_load,'\DHI_confinement_12.png']);



fnt_tit = 12;
fig31 = figure(31);
subplot(1,3,1); hold on;
imagesc(x_twin,y_twin,dhi_inverted);
h32 = scatter(x_twin,y_twin(centroid_abs),'k.');
set(h32,'sizedata',1000);
sct1 = scatter(x_twin,a_l,'m.');
set(sct1,'sizedata',1000);
sct2 = scatter(x_twin,a_r,'m.');
set(sct2,'sizedata',1000);
set(gca,'ydir','normal');
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .083]);
title('n_e [m^{-3}]','fontsize',fnt_tit);

subplot(1,3,2); hold on;
imagesc(x_twin,y_twin,B_mat);
% h32 = scatter(x_twin,y_twin(centroid_abs),'k','.');
% set(h32,'sizedata',100);
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .083]);
title('B Field [T]','fontsize',fnt_tit);
% set(gca,'clim',[0 15]);

subplot(1,3,3); hold on;
imagesc(x_twin,y_twin,T_mat);
xlim([x_twin(1) x_twin(end)]);
ylim([y_twin(1) y_twin(end)]);
colorbar;
colormap jet;
% title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
xlabel('Axial distance [m]','fontsize',fnt);
ylabel('Impact parameter [m]','fontsize',fnt);
set(gca,'fontsize',fnt);
set(gca,'xtick',[.077 .083]);
title('T_e = T_i [eV]','fontsize',fnt_tit);
% set(gca,'clim',[0 5000]);





%% subplots of normalized scaling:
scatter_colors = colormap(copper);
fig40 = figure(40);

mar_size = 200;
x_line = linspace(0,4,100);
y_line = linspace(0,4,100);

%% density vs. size adjusted for linear density:
subplot(1,3,1); hold on;

% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(1)).^(-2);
% n2_a2_vec_N1 = den_char_l(i)*(a2_vec/a_char_l(i)).^(-2)*(1);
% n2_a2_vec_N2 = n1*(a2_vec/a1).^(-2)*(1/4);
% n2_a2_vec_N3 = n1*(a2_vec/a1).^(-2)*(1/8);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2)*(N_linear_density_l(i)/N_linear_density_l(1))^(1)...
        ,(den_char_l(i)/den_char_l(1)),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2),(den_char_r(i)/den_char_l(1))...
%         *(N_linear_density_r(i)/N_linear_density_l(1))^(-1),'k','s','filled');    
%     if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

plot(x_line,y_line,'k:','linewidth',2);
% p1 = patch([0.005 0.005 0.015 0.015],[1e22 1e23 1e23 1e22],'r','EdgeColor','r');
% alpha(p1,0.25);

% plt1 = plot(a2_vec,n2_a2_vec,'k');
% plt2 = plot(a2_vec,n2_a2_vec_N1,'r');
% plt3 = plot(a2_vec,n2_a2_vec_N2,'b');
% plt4 = plot(a2_vec,n2_a2_vec_N3,'m');

xlabel('(a_2/a_1)^{-2}(N_2/N_1)','fontsize',fnt);
ylabel('(n_2/n_1) ','fontsize',fnt);
title('Density','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
lim_max = 1.6;
lim_min = .75;
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% leg1 = legend([plt1 plt2 plt3 plt4],'N_2/N_1 = 1','N_2/N_1 = 1/2','N_2/N_1 = 1/4','N_2/N_1 = 1/8');


gamma = 5/3;
%% field vs. size:
subplot(1,3,2); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
% B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-gamma)*(N_linear_density_l(i)/N_linear_density_l(1))^(gamma/2),...
        (B_char_l(i)/B_char_l(1)),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-gamma),(B_char_r(i)/B_char_l(1))*...
%         (N_linear_density_r(i)/N_linear_density_l(1))^(-exponent_1),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end

% p1 = patch([0.005 0.005 0.015 0.015],[.3 1.8 1.8 0.3],'r','EdgeColor','r');
% alpha(p1,0.25);
% % plot(a2_vec,B2_a2_vec,'k');
% plot(a2_vec,B2_a2_vec_N1,'r');
% plot(a2_vec,B2_a2_vec_N2,'b');
% plot(a2_vec,B2_a2_vec_N3,'m');
plot(x_line,y_line,'k:','linewidth',2);


xlabel('(a_2/a_1)^{-\gamma}  (N_2/N_1) ^{\gamma^2/2}','fontsize',fnt);
ylabel('(B_2/B_1)','fontsize',fnt);
title('Magnetic field','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);
xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig21,[directory_load,'\DHI_confinement_21.png']);


%% temperature vs. size:
subplot(1,3,3); hold on;

% T2_a2_vec_N1 = T_char_l(i)*(a2_vec/a_char_l(i)).^(-2*(gamma-1))*(1)^(gamma-1);
% T2_a2_vec_N2 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/4)^(gamma-1);
% T2_a2_vec_N3 = T1*(a2_vec/a1).^(-2*(gamma-1))*(1/8)^(gamma-1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses
%     sct1 = scatter(a_char_l(i),T_char_l(i)*(N_linear_density_l(1)/N_linear_density_l(i))^(gamma-1),'k','d','filled');
    sct1 = scatter((a_char_l(i)/a_char_l(1))^(-2*(gamma-1))*(N_linear_density_l(i)/N_linear_density_l(1))^(-1+gamma),...
        (T_char_l(i)/T_char_l(1)),'k','d','filled');
    if i < size(scatter_colors,1)
        set(sct1,'CData',scatter_colors(i,:));
    else
        set(sct1,'CData',scatter_colors(end,:));
    end
%     set(sct1,'MarkerEdgeColor','k');
    set(sct1,'sizedata',mar_size);
%     sct2 = scatter(a_char_r(i),T_char_r(i)*(N_linear_density_l(1)/N_linear_density_r(i))^(gamma-1),'r','s','filled');
%     sct2 = scatter((a_char_r(i)/a_char_r(1))^(-2*(gamma-1)),...
%         (T_char_r(i)/T_char_l(1))*(N_linear_density_r(i)/N_linear_density_l(1))^(1-gamma),'k','s','filled');    if i < size(scatter_colors,1)
%         set(sct2,'CData',scatter_colors(i,:));
%     else
%         set(sct2,'CData',scatter_colors(end,:));
%     end
%     set(sct2,'MarkerEdgeColor','k');
end
% p2 = patch([0.005 0.005 0.015 0.015],[20 220 220 20],'r','EdgeColor','r');
% alpha(p2,0.25);
% plot(a2_vec,T2_a2_vec,'k');
% plot(a2_vec,T2_a2_vec_N1,'r');
% plot(a2_vec,T2_a2_vec_N2,'b');
% plot(a2_vec,T2_a2_vec_N3,'m');
plot(x_line,y_line,'k:','linewidth',2);

xlabel('(a_2/a_1)^{-2(\gamma-1)}(N_2/N_1)^{\gamma-1}','fontsize',fnt);
ylabel('(T_2/T_1)','fontsize',fnt);
title('Temperature','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);


xlim([lim_min lim_max]);
ylim([lim_min lim_max]);

% saveas(fig12,[directory_load,'\DHI_confinement_12.png']);



figure(50);
mar_size = 50;

subplot(1,2,1); hold on;

exponent_1 = (gamma^2-gamma)/(2*(gamma-1));
% B2_a2_vec_N1 = B_char_l(i)*(a_char_l(i)./a2_vec).^(gamma)*(1)^(exponent_1);
% B2_a2_vec_N2 = B1*(a1./a2_vec).^(-gamma)*(1/4)^(exponent_1);
% B2_a2_vec_N3 = B1*(a1./a2_vec).^(-gamma)*(1/8)^(exponent_1);

for i = length(x_twin):-1:1 % iterate through all plasma pulses 
    ratio_B(i) = (B_char_l(i)/B_char_l(1))/(a_char_l(i)/a_char_l(1))^(-gamma)*(N_linear_density_l(i)/N_linear_density_l(1))^(gamma/2);
end

sct1 = scatter(x_twin,ratio_B,'k','filled');
set(sct1,'sizedata',mar_size);

line([min(x_twin) max(x_twin)],[1 1],'color','k','linestyle',':','linewidth',2);
xlim([min(x_twin) max(x_twin)]);
ylim([0 2]);

xlabel('Axial distance [m]','fontsize',fnt);
ylabel('(B_2/B_1)(a_2/a_1)^{\gamma}  (N_2/N_1) ^{-\gamma^2/2}','fontsize',fnt);
title('Magnetic field','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);

subplot(1,2,2); hold on;

for i = length(x_twin):-1:1 % iterate through all plasma pulses
    ratio_T(i) = (T_char_l(i)/T_char_l(1))/(a_char_l(i)/a_char_l(1))^(-2*(gamma-1))*(N_linear_density_l(i)/N_linear_density_l(1))^(-1+gamma);
end

sct2 = scatter(x_twin,ratio_T,'k','filled');
set(sct2,'sizedata',mar_size);

line([min(x_twin) max(x_twin)],[1 1],'color','k','linestyle',':','linewidth',2);
xlim([min(x_twin) max(x_twin)]);
ylim([0 2]);

xlabel('Axial distance [m]','fontsize',fnt);
ylabel('(T_2/T_1)(a_2/a_1)^{2(\gamma-1)}(N_2/N_1)^{1-\gamma}','fontsize',fnt);
title('Temperature','fontsize',fnt);
% legend([p1, p2],'ZaP','ZaP-HD');
set(gca,'fontsize',fnt);


% xlim([lim_min lim_max]);
% ylim([lim_min lim_max]);


return;


% [B_prof,T_prof,v_drift,den_char,B_char,...
%     T_char,a_char] = DHI_profiles(rad_r{k},den_num_axial_r_plt{k}',I,0.02);
