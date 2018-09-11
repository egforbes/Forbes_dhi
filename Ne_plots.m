clear all; clc;
close all;

load('Ne_plots_baseline');

x_base = x_twin_red_adj;
y_base = y_twin_red_adj;
den_base = den_int_full;
cross_base = cross_sect;

load('Ne_plots_plasma');

x_plasma = x_twin_red_adj;
y_plasma = y_twin_red_adj;
den_plasma = den_int_full;
cross_plasma = cross_sect;


denmax = 1.5e21;
fnt1 = 20;
fnt2 = 16;

figure(1); 
subplot(2,2,1); hold on;
h1 = pcolor(x_base,y_base,den_base);
set(h1,'edgecolor','none');
colormap jet
colorbar;
set(gca,'clim',[0 denmax]);
xlim([x_base(1) x_base(end)]);
ylim([y_base(1) y_base(end)]);
title('Vacuum','fontsize',fnt2);
xlabel('Axial distance [m]','fontsize',fnt2);
ylabel('Impact parameter [m]','fontsize',fnt2);
set(gca,'fontsize',fnt2);
set(gca,'xtick',[.077 .080 .083]);
line([x_base(cross_base(1)) x_base(cross_base(1))],[y_base(1) y_base(end)],'color','k','linestyle',':','linewidth',3);
if length(cross_base)==3
    line([x_base(cross_base(2)) x_base(cross_base(2))],[y_base(1) y_base(end)],'color','k','linestyle','--','linewidth',3);%,'linestyle',':');
    line([x_base(cross_base(3)) x_base(cross_base(3))],[y_base(1) y_base(end)],'color','k','linestyle','-','linewidth',3);%,'linestyle','--');
else
end

subplot(2,2,2); hold on;
h2 = pcolor(x_plasma,y_plasma,den_plasma);
set(h2,'edgecolor','none');
colormap jet
colorbar;
set(gca,'clim',[0 denmax]);
xlim([x_plasma(1) x_plasma(end)]);
ylim([y_plasma(1) y_plasma(end)]);
title('Plasma','fontsize',fnt2);
xlabel('Axial distance [m]','fontsize',fnt2);
ylabel('Impact parameter [m]','fontsize',fnt2);
set(gca,'fontsize',fnt2);
set(gca,'xtick',[.077 .080 .083]);
line([x_base(cross_base(1)) x_base(cross_base(1))],[y_base(1) y_base(end)],'color','k','linestyle',':','linewidth',3);
if length(cross_base)==3
    line([x_base(cross_base(2)) x_base(cross_base(2))],[y_base(1) y_base(end)],'color','k','linestyle','--','linewidth',3);%,'linestyle',':');
    line([x_base(cross_base(3)) x_base(cross_base(3))],[y_base(1) y_base(end)],'color','k','linestyle','-','linewidth',3);%,'linestyle','--');
else
end

subplot(2,2,3); hold on;
plot(y_base,den_base(:,round(0.25*size(den_int_full,2))),'k:','linewidth',3);
plot(y_base,den_base(:,round(0.5*size(den_int_full,2))),'k--','linewidth',3);
plot(y_base,den_base(:,round(0.75*size(den_int_full,2))),'k-','linewidth',3);
ylim([0 denmax]);
xlim([y_base(1) y_base(end)]);
xlabel('Impact parameter [m]','fontsize',fnt2);
ylabel('N_e [m^{-2}]','fontsize',fnt2);
set(gca,'fontsize',fnt2);

subplot(2,2,4); hold on;
plot(y_plasma,den_plasma(:,round(0.25*size(den_int_full,2))),'k:','linewidth',3);
plot(y_plasma,den_plasma(:,round(0.5*size(den_int_full,2))),'k--','linewidth',3);
plot(y_plasma,den_plasma(:,round(0.75*size(den_int_full,2))),'k-','linewidth',3);
ylim([0 denmax]);
xlim([y_base(1) y_base(end)]);
xlabel('Impact parameter [m]','fontsize',fnt2);
ylabel('N_e [m^{-2}]','fontsize',fnt2);
set(gca,'fontsize',fnt2);

tit1 = suptitle('Line-integrated electron density, N_e [m^{-2}]');
set(tit1,'fontsize',fnt1);





figure(2);
subplot(1,2,2); hold on;
h1 = pcolor(x_base,y_base,den_base); axis square;
set(h1,'edgecolor','none');
colormap jet
colorbar;
set(gca,'clim',[0 denmax]);
xlim([x_base(1) x_base(end)]);
ylim([y_base(1) y_base(end)]);
% title('Vacuum','fontsize',fnt2);
xlabel([{'Axial distance [m]'},{'(b)'}],'fontsize',fnt2);
ylabel('Impact parameter [m]','fontsize',fnt2);
set(gca,'fontsize',fnt2);
set(gca,'xtick',[.077 .080 .083]);


subplot(1,2,1); hold on;
h2 = pcolor(x_plasma,y_plasma,den_plasma); axis square
set(h2,'edgecolor','none');
sct1 = scatter(x_plasma,y_plasma(centroid_abs),'k','.');
set(sct1,'sizedata',400);
colormap jet
colorbar;
set(gca,'clim',[0 denmax]);
xlim([x_plasma(1) x_plasma(end)]);
ylim([y_plasma(1) y_plasma(end)]);
% title('Plasma','fontsize',fnt2);
xlabel([{'Axial distance [m]'},{'(a)'}],'fontsize',fnt2);
ylabel('Impact parameter [m]','fontsize',fnt2);
set(gca,'fontsize',fnt2);
set(gca,'xtick',[.077 .080 .083]);

% tit1 = suptitle('Line-integrated electron density, N_e [m^{-2}]');
set(tit1,'fontsize',fnt1);