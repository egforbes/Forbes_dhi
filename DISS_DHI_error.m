% DISS_DHI_error.m
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code performs the holographic reconstruction, abel
% inverts the data, and conducts error analysis of the resulting density
% profiles.  

clear all; clc;
close all;
directory_save = 'C:\Dropbox\MPROSS\DHIpaper\Coding_rebuttal';

%% User inputs, file directories, and constants:

% Importing the desired shot numbers:
% [shots,d,xmin,xmax,ymin,ymax,sign_twin] = shotnumber_inputs(2,1);

        shots = 180215000 + [12 13];
%         shots_tmp=shot_base_pair(shots,whichPair);
        d=0.55;
        xmin=0.04856;  %0.03772;
        xmax=0.06344;   %0.05058;
        ymin=0.03112; %0.02528;
        ymax=0.05413; %0.03804;



d_real = d;


% Fraction along z axis for each cross-section:
cross_sect_frac = [.25 0.5 0.75];

% Axial measurement location:
z_loc = 0.08; % [m]

% Declare constants and ZaP plasma parameters:
e_0 = 8.85e-12; % Permitivity of free space
mu_0 = 1.257e-6; % Permeability of free space
k = 1.38e-23; % Boltzmann's constant in SI units
q = 1.6e-19; % electron charge [C]
c = 3e8; % speed of light [m/s]
mi = 1.67e-27; % proton mass [kg]
me = 9.12e-31; % electron mass [kg]
lambda = 532e-9; % laser wavelength [m]
R_electrode = 8*2.54/100;
% r = 0.01;
% T = 100*11600; % degrees Kelvin
% Te = T/2;
% Ti = T/2;
% ne = 1e24;%1e22; % [m^-3]
% ni = ne;
% V = 150e3; % [m/s]
% Bw = 0.1; %[T]

% Synthetic hologram parameters:
step_sm = 3.85e-6; % pixel size of synthetic holograms [m]
d_fringe = 20e-6; % fringe size on synthetic holograms [m]
theta = asind(lambda/d_fringe);
h = d*tand(theta);

%% Fresnel transform reconstruction, Abel inversion, and error analysis:
% Iterating through the desired shots:
for shot_ind = 1:size(shots,1)
    
    shotnum_base = shots(shot_ind,2);
    shotnum_def = shots(shot_ind,1);
    
    %% Fresnel transform reconsruction:
    [unwrapped,x_twin,y_twin,phase_diff,phase_x,phase_y] = holographic_reconstruction(shotnum_base,shotnum_def,xmin,xmax,ymin,ymax,d_real);
    sign_twin = -1;
    den_int_full = sign_twin*[unwrapped/((-q^2/(4*pi*c^2*me*e_0))*lambda)]';
    den_int_full = den_int_full - min(min(den_int_full));
    
    % reducing the size of x_twin and y_twin to match unwrapped (recall that
    % unwrapped is reduced in size because of boxcar smoothing:
    del_y = y_twin(2) - y_twin(1);
    del_x = x_twin(2) - x_twin(1);
    y_twin_red = [0:del_y:del_y*(size(den_int_full,2)-1)];
    x_twin_red = [0:del_x:del_x*(size(den_int_full,1)-1)];
    
    % adjusting the zero value of x_twin and y_twin so the origin is centered.
    % Can adjust here for when the hologram is not centered at 0 impact param.
    y_twin_red_adj = (y_twin_red - y_twin_red(end)/2);
    x_twin_red_adj = x_twin_red - x_twin_red(end)/2+z_loc;
    
    % identify number of cross sections:
    num_of_cross_sect = size(unwrapped,2);
    
        % Define the indices for the selected cross sections:
    cross_sect = round(cross_sect_frac*size(unwrapped,2));
    
    % identify number of centroids:
    num_of_centroids = 10; % sets how many centroids to try
    
    % identify number of n_edge values:
    num_of_n_edge = 15; % sets how many edge densities to try
    n_edge_delta = 5e21; % sets the step between attmpeted edge density
    
    %% Identify the optimal centroid and edge density for the Abel
    % inversion for each impact parameter:
    
    % iterate through cross sections:
    for k = 1:num_of_cross_sect
        
        % identify the max Ne value and its index:
        den_int = den_int_full(k,:);
        [Ne_max(k),ind_Ne_max(k)] = max(den_int);
        
        % initialize centroid_tmp to enable determination of centroid_abs (the
        % centroid locations in absolute pixel-space:
        centroid_ind_tmp = ind_Ne_max(k)-num_of_centroids;
        
        % Case statements to allow adjustments to the analysis for centroid
        % searches starting at or near the bounds of the data:
        quit2nextShot = 0; % for debugging
        % for centroids at the lower boundary:
        if centroid_ind_tmp <=0
            den_num_l{k} = NaN(size(unwrapped,1),1);%[];
            den_num_r{k} = NaN(size(unwrapped,1),1);
            centroid_abs(k) = 1;
            centroid(k) = 0;
            % for centroids at the upper boundary:
        elseif centroid_ind_tmp >= size(unwrapped,1)
            den_num_l{k} = NaN(size(unwrapped,1),1);
            den_num_r{k} = NaN(size(unwrapped,1),1);%[];
            centroid_abs(k) = size(unwrapped,1);
            centroid(k) = 0;
        else
            
            % for centroids near (within 2*num_of_centroids of) the upper
            % boundary:
            if centroid_ind_tmp >= size(unwrapped,1)-2*num_of_centroids
                centroids_iterations = 2*num_of_centroids-...
                    (centroid_ind_tmp-(size(unwrapped,1)-2*...
                    num_of_centroids));
                % for centroids near the lower boundary:
            else
                centroids_iterations = 2*num_of_centroids;
            end
            
            % iterate through centroids:
            for i = 1:2*num_of_centroids
                
                % selecting the left side of the Ne profile and flipping it:
                Ne_1d_left{k,i} = fliplr(den_int(1:centroid_ind_tmp));
                
                % selecting the right side of the Ne profile:
                Ne_1d_right{k,i} = den_int(centroid_ind_tmp+1:end);
                
                % iterate the index to prepare the code to move to the next centroid
                % value on the next loop iteration:
                centroid_ind_track(k,i) = centroid_ind_tmp;
                centroid_ind_tmp = centroid_ind_tmp + 1;
                
                % conductig the Abel inversion to compute the ne profile on each side
                % of the assumed centroid:
                ne_1d_left{k,i} = abel_invert(y_twin,[Ne_1d_left{k,i}-min(Ne_1d_left{k,i})]');
                ne_1d_right{k,i} = abel_invert(y_twin,[Ne_1d_right{k,i}-min(Ne_1d_right{k,i})]');
                
                % when the centroid does not divide the Ne profile in half, the longer
                % of the two sides must be truncated to match the length of the shorter
                % side.  The following case statements conduct this truncation:
                len_left(k,i) = length(Ne_1d_left{k,i});
                len_right(k,i) = length(Ne_1d_right{k,i});
                
                if len_left(k,i)>len_right(k,i)
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i}(1:len_right(k,i));
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i};
                elseif len_left(k,i)<len_right(k,i)
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i};
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i}(1:len_left(k,i));
                else
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i};
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i};
                end
                
                % identify the nominal n_edge_0 value as the density in the longer
                % profile at the radius of the edge of the short profile:
                if len_left(k,i)>len_right(k,i)
                    n_edge_0(k,i) = ne_1d_left_trunc{k,i}(len_right(k,i));
                elseif len_left(k,i)<len_right(k,i)
                    n_edge_0(k,i) = ne_1d_right_trunc{k,i}(len_left(k,i));
                else
                    n_edge_0(k,i) = 0;
                end
                
                % initialize n_edge_tmp to enable determination of centroid_abs (the
                % centroid locations in absolute pixel-space:
                n_edge_ind_tmp = -num_of_n_edge;
                
                % iterate through n_edge values:
                for j = 1:2*num_of_n_edge
                    % populating list of all n_edge attempted:
                    n_edge_val(j) = n_edge_ind_tmp*n_edge_delta;
                    
                    % Iterating through n_edge values for the shorter profile:
                    if len_left(k,i)>len_right(k,i)
                        ne_1d_left_compare{k,i,j} = ne_1d_left_trunc{k,i};
                        ne_1d_right_compare{k,i,j} = ne_1d_right_trunc{k,i}+n_edge_delta*n_edge_ind_tmp;
                    elseif len_left(k,i)<len_right(k,i)
                        ne_1d_left_compare{k,i,j} = ne_1d_left_trunc{k,i}+n_edge_delta*n_edge_ind_tmp;
                        ne_1d_right_compare{k,i,j} = ne_1d_right_trunc{k,i};
                    else
                        if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                            ne_1d_left_compare{k,i,j} = ne_1d_left_trunc{k,i}+n_edge_delta*n_edge_ind_tmp;
                            ne_1d_right_compare{k,i,j} = ne_1d_right_trunc{k,i};
                        else
                            ne_1d_left_compare{k,i,j} = ne_1d_left_trunc{k,i};
                            ne_1d_right_compare{k,i,j} = ne_1d_right_trunc{k,i}+n_edge_delta*n_edge_ind_tmp;
                        end
                    end
                    % iterate the index to prepare the code to move to the next
                    % n_edge value on the next loop iteration:
                    n_edge_ind_track(k,i,j) = n_edge_ind_tmp;
                    n_edge_ind_tmp = n_edge_ind_tmp + 1;
                    
                    % computing the L2 norm of the difference in the sides of the ne profile:
%                     dens_diff{k,i,j} = abs(ne_1d_left_compare{k,i,j}(1:3)-...
%                         ne_1d_right_compare{k,i,j}(1:3));
%                     norm_dens(k,i,j) = norm(dens_diff{k,i,j},2)/...
%                         (length(dens_diff{k,i,j}(1:3)));
                    dens_diff{k,i,j} = abs(ne_1d_left_compare{k,i,j} - ne_1d_right_compare{k,i,j});
                    norm_dens(k,i,j) = norm(dens_diff{k,i,j},2)/(length(dens_diff{k,i,j})); %(1:length(dens_diff{k,i,j}))
                    % norm_dens(k,i,j) = norm(dens_diff{k,i,j}/(length(dens_diff{k,i,j})-1),2); %(1:length(dens_diff{k,i,j}))
                    
                end
                
                % changing all zero-valued elements in norm_dens to NaNs (allows the code to ignore them when searching for the minimum)
                ind_zeros = find(norm_dens(k,i,:)==0);
                norm_dens(k,i,ind_zeros) = NaN;
                
                % find the minimum L2 norm value for each combination of cross-section and centroid:
                [norm_intermediate(k,i),n_edge_intermediate(k,i)] = min(norm_dens(k,i,:),[],3); % n_edge_intermediate records the index location of the n_edge yielding minimum L2 norm
                
            end
            
            % find the minimum L2 norm value for each cross-section:
            [norm_min_val(k),centroid(k)] = min(norm_intermediate(k,:),[],2); % centroid records the relative index location of the centroid location yielding minimum L2 norm
            n_edge(k) = n_edge_intermediate(k,centroid(k));
            
            % assigning the final left and right number density profiles:
            if len_left(k,centroid(k))>len_right(k,centroid(k))
                den_num_l{k} = ne_1d_left{k,centroid(k)};
                den_num_r{k} = ne_1d_right{k,centroid(k)}+n_edge_delta*(n_edge(k)-num_of_n_edge);
                den_int_edge_axial(k) = den_int_full(k,1);
                
            elseif len_left(k,centroid(k))<len_right(k,centroid(k))
                den_num_l{k} = ne_1d_left{k,centroid(k)}+n_edge_delta*(n_edge(k)-num_of_n_edge);
                den_num_r{k} = ne_1d_right{k,centroid(k)};
                den_int_edge_axial(k) = den_int_full(k,end);
                
            else
                if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                    den_num_l{k} = ne_1d_left_trunc{k,centroid(k)}+n_edge_delta*(n_edge(k)-num_of_n_edge);
                    den_num_r{k} = ne_1d_right_trunc{k,centroid(k)};
                    den_int_edge_axial(k) = den_int_full(k,end);
                    
                else
                    den_num_l{k} = ne_1d_left_trunc{k,centroid(k)};
                    den_num_r{k} = ne_1d_right_trunc{k,centroid(k)}+n_edge_delta*(n_edge(k)-num_of_n_edge);
                    den_int_edge_axial(k) = den_int_full(k,end);
                    
                end
            end
            
            % computing the centroid indices in absoluter pixel coordinates:
            centroid_abs(k) = centroid_ind_track(k,centroid(k));
            
        end
    end
        
        if quit2nextShot == 0
            
            %% 4.) Forming the den_num_full matrix:
            den_num_full = zeros(size(den_int_full));
            for k = 1:num_of_cross_sect
                den_num_full(k,(1+centroid_abs(k)-length(den_num_l{k})):centroid_abs(k)) = flipud(den_num_l{k});
                den_num_full(k,(1+centroid_abs(k)):(centroid_abs(k)+length(den_num_r{k}))) = den_num_r{k};
            end
            
            %% Computing linear density:
            for k = 1:num_of_cross_sect
                rad_l{k} = 0:del_y:(length(den_num_l{k})-1)*del_y;  % del_y in empirical reconstruction coords.
                rad_r{k} = 0:del_y:(length(den_num_r{k})-1)*del_y; % del_y in empirical reconstruction coords.
                
                N_linear_density_l_tmp = cumtrapz(rad_l{k},2*pi*den_num_l{k}.*rad_l{k}');
                N_linear_density_l(k) = N_linear_density_l_tmp(end);
                
                N_linear_density_r_tmp = cumtrapz(rad_r{k},2*pi*den_num_r{k}.*rad_r{k}');
                N_linear_density_r(k) = N_linear_density_r_tmp(end);
            end
            clear rad_l rad_r
            
            %% Preparing to plot:
            x_0 = x_twin_red_adj(1);
            y_0 = y_twin_red_adj(1);
            
            den_num_full = den_num_full';
            den_int_full = den_int_full';
            
            %% Plotting:
            fnt = 24;
            % plotting Ne:
            fig1 = figure(1);hold on;
            h11 = pcolor(x_twin_red_adj,y_twin_red_adj,den_int_full);%shading 'interp';
            set(h11,'edgecolor','none');
            %             h12 = scatter(x_twin_red_adj,y_twin_red_adj(centroid_abs),'k','.');
            %             set(h12,'sizedata',1000);
            colormap jet
            set(fig1,'position',[1290         520         560         420]);
            colorbar;
            xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
            ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
            % line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);
            % if length(cross_sect)==3
            %     line([cross_sect(2)*del_x+x_0 cross_sect(2)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle',':');
            %     line([cross_sect(3)*del_x+x_0 cross_sect(3)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle','--');
            % else
            % end
            title('Line-integrated electron density, N_e [m^{-2}]','fontsize',fnt);
            % title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)],'fontsize',fnt);
            xlabel('Axial distance [m]','fontsize',fnt);
            ylabel('Impact parameter [m]','fontsize',fnt);
            set(gca,'fontsize',fnt);
            set(gca,'xtick',[.077 .080 .083]);
            
            den_num_full = den_num_full - min(min(den_num_full));
            fig3 = figure(3);hold on;
            h31 = pcolor(x_twin_red_adj,y_twin_red_adj,den_num_full); %note I enforce no negative densities here
            set(h31,'edgecolor','none');
            h32 = scatter(x_twin_red_adj,y_twin_red_adj(centroid_abs),'k','.');
            set(h32,'sizedata',1000);
            colormap jet
            set(fig3,'position',[1867         519         560         420]);
            colorbar;
            xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
            ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
            line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k');
            line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);
            if length(cross_sect)==3
                line([cross_sect(2)*del_x+x_0 cross_sect(2)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle',':');
                line([cross_sect(3)*del_x+x_0 cross_sect(3)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle','--');
            else
            end
            title('Electron number density, n_e [m^{-3}]','fontsize',fnt);
            % title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
            xlabel('Axial distance [m]','fontsize',fnt);
            ylabel('Impact parameter [m]','fontsize',fnt);
            set(gca,'fontsize',fnt);
            set(gca,'xtick',[.077 .080 .083]);
            
            if length(cross_sect)==3
                
                fig4 = figure(4);hold on;
                subplot(3,1,1); hold on;
                % plot(den_num_full(:,cross_sect(1)));
                plot(den_num_l{cross_sect(1)},'color','k');
                plot(den_num_r{cross_sect(1)},'color','r');
                %     ylim([0 2e23]);
                ylabel('upstream ne');
                
                % scatter(centroid_abs(cross_sect(1)),den_num_full(centroid_abs(cross_sect(1)),cross_sect(1)),'k');
                subplot(3,1,2); hold on;
                % plot(den_num_full(:,cross_sect(2)));
                plot(den_num_l{cross_sect(2)},'color','k');
                plot(den_num_r{cross_sect(2)},'color','r');
                %     ylim([0 2e23]);
                ylabel('middle ne');
                
                % scatter(centroid_abs(cross_sect(2)),den_num_full(centroid_abs(cross_sect(2)),cross_sect(2)),'k');
                subplot(3,1,3); hold on;
                % plot(den_num_full(:,cross_sect(3)));
                plot(den_num_l{cross_sect(3)},'color','k');
                plot(den_num_r{cross_sect(3)},'color','r');
                %     ylim([0 2e23]);
                ylabel('downstream ne');
                % scatter(centroid_abs(cross_sect(3)),den_num_full(centroid_abs(cross_sect(3)),cross_sect(3)),'k');
                
                set(fig4,'position',[1868          10         560         420]);
                
            else
            end
            
            %% Conducting error analysis on three selected cross-sections:
            
            % conduct the error analysis only if both side lengths are
            % sufficiently long:
            
            for k = 1:length(cross_sect)
                side_len_check(k,:) = [length(den_num_l{cross_sect(k)}) length(den_num_r{cross_sect(k)})];
            end
            
            if min(min(side_len_check)) > 0.1*size(den_int_full,1)
                
                % setting up the radial component vector for the three cross-sections for error analysis:
                for k = 1:length(cross_sect)
                    
                    
                    n_prof_tmp_l{k} = den_num_l{cross_sect(k)};
                    n_prof_tmp_r{k} = den_num_r{cross_sect(k)};
                    
                    rad_l{k} = 0:del_y:(length(n_prof_tmp_l{k})-1)*del_y;  % del_y in empirical reconstruction coords.
                    rad_r{k} = 0:del_y:(length(n_prof_tmp_r{k})-1)*del_y; % del_y in empirical reconstruction coords.
                    
                    % need to extract from den_num_full to accout for minimum in 2d density map:
                    n_prof_l{k} = den_num_l{cross_sect(k)};
                    n_prof_r{k} = den_num_r{cross_sect(k)};
                    
                    
                    if len_left(cross_sect(k),centroid(cross_sect(k)))>len_right(cross_sect(k),centroid(cross_sect(k)))%% replaced k with cross_sect(k)
                        %         n_prof_l{k} = n_prof_l{k};
                        %         n_prof_r{k} = n_prof_r{k};%+n_edge_val(n_edge(cross_sect(k)));
                        
                        % Adding N_edge on after abel transform:
                        del_r = rad_l{k}(2)-rad_l{k}(1);
                        r_pix = len_left(cross_sect(k),centroid(cross_sect(k)));
                        y_pix = len_right(cross_sect(k),centroid(cross_sect(k)));
                        x(k) = sqrt(r_pix.^2-y_pix.^2)*del_r;
                        N_additional(k) = 2*n_prof_l{k}(end)*x(k);
                        N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}'-min(n_prof_l{k}),0)+N_additional(k);
                        N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}'-min(n_prof_r{k}),0)+N_additional(k); % don't need N_additional added here?
                        
                        % Adjusting the size of the Abel inversion matrix:
                        %         abel_diff(k) = len_left(cross_sect(k),centroid(cross_sect(k)))-len_right(cross_sect(k),centroid(cross_sect(k)));%%
                        %         N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}',0);
                        %         N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}',abel_diff(k));
                        
                    elseif len_left(cross_sect(k),centroid(cross_sect(k)))<len_right(cross_sect(k),centroid(cross_sect(k)))%%
                        %         n_prof_l{k} = n_prof_l{k};%+n_edge_val(n_edge(cross_sect(k)));
                        %         n_prof_r{k} = n_prof_r{k};
                        
                        % Adding N_edge on after:
                        del_r = rad_l{k}(2)-rad_l{k}(1);
                        r_pix = len_right(cross_sect(k),centroid(cross_sect(k)));
                        y_pix = len_left(cross_sect(k),centroid(cross_sect(k)));
                        x(k) = 2*sqrt(r_pix.^2-y_pix.^2)*del_r;
                        N_additional(k) = n_prof_r{k}(end)*x(k);
                        N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}'-min(n_prof_l{k}),0)+N_additional(k); % don't need N_additional added here?
                        N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}'-min(n_prof_r{k}),0)+N_additional(k);
                        
                        % Adjusting the size of the Abel inversion matrix:
                        %         abel_diff(k) = len_right(cross_sect(k),centroid(cross_sect(k)))-len_left(cross_sect(k),centroid(cross_sect(k)));%%
                        %         N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}',abel_diff(k));
                        %         N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}',0);
                        
                    else
                        if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                            del_r = rad_l{k}(2)-rad_l{k}(1);
                            r_pix = len_left(cross_sect(k),centroid(cross_sect(k)));
                            y_pix = len_right(cross_sect(k),centroid(cross_sect(k)));
                            x(k) = sqrt(r_pix.^2-y_pix.^2)*del_r;
                            N_additional(k) = 2*n_prof_r{k}(end)*x(k);
                            N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}'-min(n_prof_l{k}),0)+N_additional(k);
                            N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}'-min(n_prof_r{k}),0)+N_additional(k); % don't need N_additional added here?
                            
                            
                        else
                            del_r = rad_l{k}(2)-rad_l{k}(1);
                            r_pix = len_left(cross_sect(k),centroid(cross_sect(k)));
                            y_pix = len_right(cross_sect(k),centroid(cross_sect(k)));
                            x(k) = sqrt(r_pix.^2-y_pix.^2)*del_r;
                            N_additional(k) = 2*n_prof_l{k}(end)*x(k);
                            N_prof_l{k}(1:length(n_prof_l{k}),1) = abel_red(rad_l{k},n_prof_l{k}'-min(n_prof_l{k}),0)+N_additional(k);
                            N_prof_r{k}(1:length(n_prof_r{k}),1) = abel_red(rad_r{k},n_prof_r{k}'-min(n_prof_r{k}),0)+N_additional(k); % don't need N_additional added here?
                            
                            
                        end
                        
                    end
                    
                    N_edge_l(k) = N_prof_l{k}(end);
                    N_edge_r(k) = N_prof_r{k}(end);
                    %     N_prof_l{k}(1:length(n_prof_l{k}),1) = abel(rad_l{k},n_prof_l{k}');
                    %     N_prof_r{k}(1:length(n_prof_r{k}),1) = abel(rad_r{k},n_prof_r{k}');
                    
                end
                
                % conducting error analysis for each cross-section
                for k = 1:length(cross_sect)
                    
                    % converting the Ne to phase shift for each side of profile:
                    ph_shift_l_short{k} = (-q^2/(4*pi*c^2*me*e_0))*lambda*N_prof_l{k};
                    ph_shift_r_short{k} = (-q^2/(4*pi*c^2*me*e_0))*lambda*N_prof_r{k};
                    %     ph_shift_l_short{k} = ph_shift_l_short{k}-min(ph_shift_l_short{k});
                    %     ph_shift_r_short{k} = ph_shift_r_short{k}-min(ph_shift_r_short{k});
                    
                    % increase the number of pts in the profile to provide resolution in
                    % holograms to enable holographic reconstruction:
                    rad_hol_l{k} = [0:step_sm:rad_l{k}(end)];%linspace(0,rad_l{i}(end),2000);
                    rad_hol_r{k} = [0:step_sm:rad_r{k}(end)];%linspace(0,rad_r{i}(end),2000);
                    
                    % recall each N_prof is forced to have the same length with zero padding:
                    ph_shift_l{k} = interp1(rad_l{k},ph_shift_l_short{k},rad_hol_l{k},'spline');
                    ph_shift_r{k} = interp1(rad_r{k},ph_shift_r_short{k},rad_hol_r{k},'spline');
                    
                    % Applying the phase shift caused by the assumed density profile:
                    
                    % identify lengths for population of synthetic data:
                    hol_size_l = length(ph_shift_l{k});
                    hol_size_r = length(ph_shift_r{k});
                    
                    % arranging the phase into a 2-D array for each side of profile:
                    ph_2d_l = ones(hol_size_l,hol_size_l);
                    for jj1 = 1:hol_size_l
                        ph_2d_l(jj1,:) = ph_shift_l{k};
                    end
                    
                    ph_2d_r = ones(hol_size_r,hol_size_r);
                    for jj2 = 1:hol_size_r
                        ph_2d_r(jj2,:) = ph_shift_r{k};
                    end
                    
                    % Applying the phase shift caused by scene/reference misalignment:
                    linear_phase_sh = 2*pi/d_fringe;
                    
                    % setting up the 1-D phase distribution:
                    ph_lin_1d_l = linspace(0,linear_phase_sh*step_sm*hol_size_l,hol_size_l);
                    ph_lin_1d_r = linspace(0,linear_phase_sh*step_sm*hol_size_r,hol_size_r);
                    
                    % putting the 1-D phase into a 2-D array for each side of profile:
                    ph_lin_2d_l = ones(hol_size_l,hol_size_l);
                    ph_lin_2d_l(:,1) = ph_lin_1d_l';
                    for jj11 = 2:hol_size_l
                        ph_lin_2d_l(:,jj11) = ph_lin_1d_l;
                    end
                    
                    ph_lin_2d_r = ones(hol_size_r,hol_size_r);
                    ph_lin_2d_r(:,1) = ph_lin_1d_r';
                    for jj21 = 2:hol_size_r
                        ph_lin_2d_r(:,jj21) = ph_lin_1d_r;
                    end
                    
                    % summing the phase contributions into a 2-D array for each side of profile:
                    ph_tot_l = ph_2d_l + ph_lin_2d_l + .02*rand(size(ph_2d_l));% + ph_lin_lg_2d;
                    ph_base_l = ph_lin_2d_l + .02*rand(size(ph_2d_l));
                    ph_tot_r = ph_2d_r + ph_lin_2d_r + .02*rand(size(ph_2d_r));% + ph_lin_lg_2d;
                    ph_base_r = ph_lin_2d_r + .02*rand(size(ph_2d_r));
                    
                    % defining the deformed/baseline holograms for each side of profile:
                    hol_def_l{k} = cos(ph_tot_l);
                    hol_base_l{k}  = cos(ph_base_l);
                    
                    hol_def_r{k} = cos(ph_tot_r);
                    hol_base_r{k}  = cos(ph_base_r);
                    
                    % setting up hologram coordinates:
                    L_l = size(hol_def_l{k});
                    M_l = L_l(2);  % x-direction
                    N_l = L_l(1);  % y-direction
                    x_hol_l{k} = 0:step_sm:(M_l-1)*step_sm;
                    y_hol_l{k} = 0:step_sm:(N_l-1)*step_sm;
                    
                    L_r = size(hol_def_r{k});
                    M_r = L_r(2);  % x-direction
                    N_r = L_r(1);  % y-direction
                    x_hol_r{k} = 0:step_sm:(M_r-1)*step_sm;
                    y_hol_r{k} = 0:step_sm:(N_r-1)*step_sm;
                    
                    % making the holograms match the format of the jpegs output from camera:
                    
                    % enforce no non-zero values of intensity:
                    hol_base_1_l = hol_base_l{k} - min(min(hol_base_l{k}));
                    hol_def_1_l = hol_def_l{k} - min(min(hol_def_l{k}));
                    
                    hol_base_1_r = hol_base_r{k} - min(min(hol_base_r{k}));
                    hol_def_1_r = hol_def_r{k} - min(min(hol_def_r{k}));
                    
                    % set the holograms to have an 8-bit depth:
                    hol_base_l{k} = (255/max(max(hol_base_1_l)))*hol_base_1_l;
                    hol_def_l{k} = (255/max(max(hol_def_1_l)))*hol_def_1_l;
                    
                    hol_base_r{k} = (255/max(max(hol_base_1_r)))*hol_base_1_r;
                    hol_def_r{k} = (255/max(max(hol_def_1_r)))*hol_def_1_r;
                    
                    hol_base_l{k}(hol_base_l{k}>255) = 255;
                    hol_base_l{k}(hol_base_l{k}<0) = 0;
                    hol_def_l{k}(hol_def_r{k}>255) = 255;
                    hol_def_l{k}(hol_def_r{k}<0) = 0;
                    
                    hol_base_r{k}(hol_base_l{k}>255) = 255;
                    hol_base_r{k}(hol_base_l{k}<0) = 0;
                    hol_def_r{k}(hol_def_r{k}>255) = 255;
                    hol_def_r{k}(hol_def_r{k}<0) = 0;
                    
                    % conduct the holographic reconstruction for each side of profile:
                    d = .5;
                    [phase_diff_l{k},x_l{k},y_l{k},xsize_fres_l(k),ysize_fres_l(k)] = DHI_rec(hol_base_l{k},hol_def_l{k},d,step_sm);
                    [phase_diff_r{k},x_r{k},y_r{k},xsize_fres_r(k),ysize_fres_r(k)] = DHI_rec(hol_base_r{k},hol_def_r{k},d,step_sm);
                    
                    % extracting the twin image:
                    
                    % set up the coordinates of the reconstruction plane:
                    x_plt_l{k} = -max(x_l{k})/2:x_l{k}(2)-x_l{k}(1):max(x_l{k})/2;
                    y_plt_l{k} = -max(y_l{k})/2:y_l{k}(2)-y_l{k}(1):max(y_l{k})/2;
                    
                    x_plt_r{k} = -max(x_r{k})/2:x_r{k}(2)-x_r{k}(1):max(x_r{k})/2;
                    y_plt_r{k} = -max(y_r{k})/2:y_r{k}(2)-y_r{k}(1):max(y_r{k})/2;
                    
                    % user-defined boundary of twin image:
                    %     xmin = .07913;
                    %     xmax = .09043;
                    %     ymin = .09414;
                    %     ymax = .1053;
                    
                    [x_twin_l{k},y_twin_l{k},twin_img_l{k}] = twin_img_extract(x_plt_l{k},y_plt_l{k},phase_diff_l{k},x_hol_l{k},y_hol_l{k},h);
                    [x_twin_r{k},y_twin_r{k},twin_img_r{k}] = twin_img_extract(x_plt_r{k},y_plt_r{k},phase_diff_r{k},x_hol_r{k},y_hol_r{k},h);
                    %     return;
                    % smoothing and unwrapping the twin image:
                    width = 4; % set the width of the boxcar filter for smoothing
                    [unwrapped_l{k},phase_f,unsmoothed_unwrapped_l{k}] = smooth_unwrap_1d(twin_img_l{k},width);
                    [unwrapped_r{k},phase_f,unsmoothed_unwrapped_r{k}] = smooth_unwrap_1d(twin_img_r{k},width);
                    
                    % converting the unwrapped phase to profiles of Ne and ne (different cases assign abel_diff correctly):
                    [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,N_additional(k),x(k));
                    [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,N_additional(k),x(k));
                    %                     if len_left(k,centroid(cross_sect(k)))>len_right(k,centroid(cross_sect(k)))%%
                    %                         [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,0,0);
                    %                         [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,N_additional(k),x(k));
                    %
                    %                         %             [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,N_edge_l(k));
                    %                         %             [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),abel_diff(k),N_edge_r(k));
                    %
                    %                     elseif len_left(k,centroid(cross_sect(k)))<len_right(k,centroid(cross_sect(k)))%%
                    %                         [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,N_additional(k),x(k));
                    %                         [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,0,0);
                    %
                    %                         %             [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),abel_diff(k),N_edge_l(k));
                    %                         %             [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,N_edge_r(k));
                    %
                    %                     else
                    %
                    %                         if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                    %                             [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,N_additional(k),x(k));
                    %                             [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,0,0);
                    %
                    %                         else
                    %                             [rad_rec_l{k},den_int_rec_l{k},den_num_rec_l{k},ph_rec_l{k},A_l] = phase2profiles(unwrapped_l{k},ysize_fres_l(k),ysize_fres_l(k),0,0,0);
                    %                             [rad_rec_r{k},den_int_rec_r{k},den_num_rec_r{k},ph_rec_r{k},A_r] = phase2profiles(unwrapped_r{k},xsize_fres_r(k),ysize_fres_r(k),0,N_additional(k),x(k));
                    %
                    %                         end
                    %
                    %                     end
                end
                
                
                
                %% 6.) interpolating synthetically reconsructed density and computing error:
                for k = 1:length(cross_sect)
                    den_num_rec_l_plt{k} = interp1(rad_rec_l{k},den_num_rec_l{k},rad_l{k});
                    den_num_rec_r_plt{k} = interp1(rad_rec_r{k},den_num_rec_r{k},rad_r{k});
                    
                    n_error_l{k} = n_prof_l{k}-den_num_rec_l_plt{k}';
                    n_error_r{k} = n_prof_r{k}-den_num_rec_r_plt{k}'; %%den_numrec
                end
                
                
                
                
                
                
                
                
                
                
                
                
                
                fig401 = figure(401);
                
                k = 1;
                subplot(2,1,1); hold on; grid on;
                h1 = plot(n_prof_l{k},'color','k');
                h2 = plot(den_num_rec_l_plt{k}','color','m');
                title('left');
                subplot(2,1,2); hold on; grid on;
                plot(n_prof_r{k},'color','k');
                plot(den_num_rec_r_plt{k}','color','m');
                title('right');
                legend([h1 h2],'measured profile','after synth recon');
                
                %         figure(402);
                %         subplot(2,1,1); hold on; grid on;
                %         h1 = plot(N_prof_l{k},'color','k');
                %         h2 = plot(den_int_rec_l{k},'color','m');
                %
                %         subplot(2,1,2); hold on; grid on;
                %         plot(N_prof_r{k},'color','k');
                %         plot(den_int_rec_r{k},'color','m');
                %         legend([h1 h2],'measured profile','after synth recon');
                
                
                fig400 = figure(400);
                fnt = 18;
                fnt_2 = 14;
                fnt_3 = 18;
                col_sec = 'r';
                subplot(3,1,1);hold on;
                p1 = plot(rad_l{1},n_prof_l{1},col_sec,'linewidth',4);
                p2 = plot(rad_r{1},n_prof_r{1},'k','linewidth',4);
                % shadedErrorBar(rad_r{1},n_prof_r{1},n_error_r{1},col_sec,1);
                % shadedErrorBar(rad_r{1},n_prof_r{1},n_error_r{1},'k',1);
                h40011 = errorbar(rad_l{1},n_prof_l{1},n_error_l{1},'color',col_sec);
                h40012 = errorbar(rad_r{1},n_prof_r{1},n_error_r{1},'color','k');
                
                % p1_syn = plot(rad_l{1},den_num_rec_l_plt{1},'k','linewidth',2,'linestyle',':');
                % p2_syn = plot(rad_r{1},den_num_rec_r_plt{1},'b','linewidth',2,'linestyle',':');
                
                %                 ylim([-5e21 max([max(n_prof_l{1}),max(n_prof_r{1}),max(den_num_rec_l{1}),max(den_num_rec_r{1}) max(n_prof_l{2}),max(n_prof_r{2}),max(den_num_rec_l{2}),max(den_num_rec_r{2}) max(n_prof_l{3}),max(n_prof_r{3}),max(den_num_rec_l{3}),max(den_num_rec_r{3})])]);
                %                 xlim([rad_rec_l{1}(1) max([rad_rec_l{1}(end) rad_rec_l{2}(end) rad_rec_l{3}(end) rad_rec_r{1}(end) rad_rec_r{2}(end) rad_rec_r{3}(end)])]);
                set(gca,'fontsize',fnt_3);
                title(['z = 0.145 m'],'fontsize',fnt_2);
                % title(['Z = ',num2str(x_twin_red_adj(cross_sect(1))),' m'],'fontsize',fnt_2);
                ylabel('n_e [m^{-3}]','fontsize',fnt);
                % title(['norm = ',num2str(norm_dens(cross_sect(1),centroid(cross_sect(1))))]);
                % ylabel(['z = ', num2str(x_twin_red_adj(cross_sect(1)))]);
                % legend([p1 p2],'left','right','fontsize',fnt);%'synth left','synth right');
                
                if length(cross_sect) ==3
                    subplot(3,1,2);hold on;
                    plot(rad_l{2},n_prof_l{2},col_sec,'linewidth',4);
                    plot(rad_r{2},n_prof_r{2},'k','linewidth',4);
                    %     shadedErrorBar(rad_r{2},n_prof_r{2},n_error_r{2},col_sec,1);
                    %     shadedErrorBar(rad_r{2},n_prof_r{2},n_error_r{2},'k',1);
                    h40021 = errorbar(rad_l{2},n_prof_l{2},n_error_l{2},'color',col_sec);
                    errorbar(rad_r{2},n_prof_r{2},n_error_r{2},'color','k');
                    
                    % plot(rad_l{1},den_num_rec_l_plt{2},'k','linewidth',2,'linestyle',':');
                    % plot(rad_r{1},den_num_rec_r_plt{2},'b','linewidth',2,'linestyle',':');
                    
                    ylim([0 max([max(n_prof_l{1}),max(n_prof_r{1}),max(den_num_rec_l{1}),max(den_num_rec_r{1}) max(n_prof_l{2}),max(n_prof_r{2}),max(den_num_rec_l{2}),max(den_num_rec_r{2}) max(n_prof_l{3}),max(n_prof_r{3}),max(den_num_rec_l{3}),max(den_num_rec_r{3})])]);
                    xlim([rad_rec_l{1}(1) max([rad_rec_l{1}(end) rad_rec_l{2}(end) rad_rec_l{3}(end) rad_rec_r{1}(end) rad_rec_r{2}(end) rad_rec_r{3}(end)])]);
                    set(gca,'fontsize',fnt_3);
                    title(['z = 0.15 m'],'fontsize',fnt_2);
                    %     title(['Z = ',num2str(x_twin_red_adj(cross_sect(2))),' m'],'fontsize',fnt_2);
                    ylabel('n_e [m^{-3}]','fontsize',fnt);
                    %     title(['norm = ',num2str(norm_dens(cross_sect(2),centroid(cross_sect(2))))]);
                    %     ylabel(['z = ', num2str(x_twin_red_adj(cross_sect(2)))]);
                    
                    subplot(3,1,3);hold on;
                    plot(rad_l{3},n_prof_l{3},col_sec,'linewidth',4);
                    plot(rad_r{3},n_prof_r{3},'k','linewidth',4);
                    %     shadedErrorBar(rad_l{3},n_prof_l{3},n_error_l{3},col_sec,1);
                    %     shadedErrorBar(rad_r{3},n_prof_r{3},n_error_r{3},'k',1);
                    errorbar(rad_l{3},n_prof_l{3},n_error_l{3},'color',col_sec);
                    errorbar(rad_r{3},n_prof_r{3},n_error_r{3},'color','k');
                    
                    % plot(rad_l{1},den_num_rec_l_plt{3},'k','linewidth',2,'linestyle',':');
                    % plot(rad_r{1},den_num_rec_r_plt{3},'b','linewidth',2,'linestyle',':');
                    
                    ylim([0 max([max(n_prof_l{1}),max(n_prof_r{1}),max(den_num_rec_l{1}),max(den_num_rec_r{1}) max(n_prof_l{2}),max(n_prof_r{2}),max(den_num_rec_l{2}),max(den_num_rec_r{2}) max(n_prof_l{3}),max(n_prof_r{3}),max(den_num_rec_l{3}),max(den_num_rec_r{3})])]);
                    xlim([rad_rec_l{1}(1) max([rad_rec_l{1}(end) rad_rec_l{2}(end) rad_rec_l{3}(end) rad_rec_r{1}(end) rad_rec_r{2}(end) rad_rec_r{3}(end)])]);
                    set(gca,'fontsize',fnt_3);
                    title(['z = 0.155 m'],'fontsize',fnt_2);
                    %     title(['Z = ',num2str(x_twin_red_adj(cross_sect(3))),' m'],'fontsize',fnt_2);
                    ylabel('n_e [m^{-3}]','fontsize',fnt);
                    xlabel('radius [m]','fontsize',fnt);
                    %     title(['norm = ',num2str(norm_dens(cross_sect(3),centroid(cross_sect(3))))]);
                    %     ylabel(['z = ', num2str(x_twin_red_adj(cross_sect(3)))]);
                    
                else
                end
                
                
                if sign_twin > 0
                    %                     saveas(fig400,[directory_save,'\ne_crosssect_',num2str(shotnum_def),'_',num2str(shotnum_base),'POS.tif']);
                    %                     save([directory_save,'\ne_cross_sections_',num2str(shotnum_def),'_',num2str(shotnum_base),'.mat'],'sign_twin','rad_l','rad_r','n_prof_l','n_prof_r','n_error_l','n_error_r');
                else
                    %                     saveas(fig400,[directory_save,'\ne_crosssect_',num2str(shotnum_def),'_',num2str(shotnum_base),'.tif']);
                    %                     save([directory_save,'\ne_cross_sections_',num2str(shotnum_def),'_',num2str(shotnum_base),'.mat'],'sign_twin','rad_l','rad_r','n_prof_l','n_prof_r','n_error_l','n_error_r');
                end
                
            else % this else is for choosing whether or not to try the error analysis.
            end % this end is for choosing whether or not to try the error analysis.
            
            
            %% 7.) accounting for axial variation:
            
            for k = 1:num_of_cross_sect
                y_edge_l{k} = 0:del_y:(length(den_num_l{k})-1)*del_y;  % del_y in empirical reconstruction coords.
                y_edge_r{k} = 0:del_y:(length(den_num_r{k})-1)*del_y;
                
                y_count(k) = max(y_edge_l{k});
                x_edge_l(k) = 2*sqrt(R_electrode^2-max(y_edge_l{k})^2);
                x_edge_r(k) = 2*sqrt(R_electrode^2-max(y_edge_r{k})^2);
                
                if len_left(k,centroid(k))>len_right(k,centroid(k))%%
                    den_num_edge_axial_l(k) = den_int_edge_axial(k)/x_edge_l(k);
                    den_num_edge_axial_r(k) = 0;
                    
                    den_num_axial_l{k} = den_num_l{k} + den_num_edge_axial_l(k);
                    den_num_axial_r{k} = den_num_r{k} + den_num_edge_axial_l(k);
                    
                elseif len_left(k,centroid(k))<len_right(k,centroid(k))%%
                    den_num_edge_axial_l(k) = 0;
                    den_num_edge_axial_r(k) = den_int_edge_axial(k)/x_edge_r(k);
                    
                    den_num_axial_l{k} = den_num_l{k} + den_num_edge_axial_r(k);
                    den_num_axial_r{k} = den_num_r{k} + den_num_edge_axial_r(k);
                    
                else
                    
                    if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                        den_num_edge_axial_l(k) = den_int_edge_axial(k)/x_edge_l(k);
                        den_num_edge_axial_r(k) = 0;
                        
                        den_num_axial_l{k} = den_num_l{k} + den_num_edge_axial_l(k);
                        den_num_axial_r{k} = den_num_r{k} + den_num_edge_axial_l(k);
                        
                    else
                        den_num_edge_axial_l(k) = 0;
                        den_num_edge_axial_r(k) = den_int_edge_axial(k)/x_edge_r(k);
                        
                        den_num_axial_l{k} = den_num_l{k} + den_num_edge_axial_r(k);
                        den_num_axial_r{k} = den_num_r{k} + den_num_edge_axial_r(k);
                    end
                    
                    
                end
                
                %                     den_num_l{k} = den_num_l{k} + den_num_edge_axial_l(k);
                %                     den_num_r{k} = den_num_r{k} + den_num_edge_axial_r(k);
                
            end
            
            % den_num_full_axial = zeros(size(den_int_full));
            for k = 1:num_of_cross_sect
                den_num_full_axial(k,(1+centroid_abs(k)-length(den_num_axial_l{k})):centroid_abs(k)) = flipud(den_num_axial_l{k});
                den_num_full_axial(k,(1+centroid_abs(k)):(centroid_abs(k)+length(den_num_axial_r{k}))) = den_num_axial_r{k};
            end
            
            min_den = min(min(den_num_full_axial));
            den_num_full_axial = den_num_full_axial'-min_den;
            
            for k = 1:length(cross_sect)
                den_num_axial_l_plt{k} = flipud(den_num_full_axial(1:length(den_num_l{cross_sect(k)}),cross_sect(k)));
                den_num_axial_r_plt{k} = den_num_full_axial((length(den_num_l{cross_sect(k)})+1):end,cross_sect(k));
            end
            
            
            
            %%% Try saving here!!!
            save([directory_save,'\ne_cross_sections_',num2str(shotnum_def),'_',num2str(shotnum_base),'.mat'],'sign_twin','rad_l','rad_r','den_num_axial_l_plt','den_num_axial_r_plt','n_error_l','n_error_r');
            
            
            
            fig33 = figure(33);hold on; fnt = 24;
            h31 = pcolor(x_twin_red_adj,y_twin_red_adj,den_num_full_axial);%shading 'interp';
            set(h31,'edgecolor','none');
            h32 = scatter(x_twin_red_adj,y_twin_red_adj(centroid_abs),'k','.');
            set(h32,'sizedata',1000);
            colormap jet
            set(fig33,'position',[1867         519         560         420]);
            colorbar;
            xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
            ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
            %             line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k');
            line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);
            %             if length(cross_sect)==3
            %                 line([cross_sect(2)*del_x+x_0 cross_sect(2)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle',':');
            %                 line([cross_sect(3)*del_x+x_0 cross_sect(3)*del_x+x_0],[y_twin_red_adj(1) y_twin_red_adj(end)],'color','k','linewidth',3);%,'linestyle','--');
            %             else
            %             end
            title('Electron number density, n_e [m^{-3}]','fontsize',fnt);
            % title(['def = ',num2str(shotnum_def),'; base = ',num2str(shotnum_base)]);
            xlabel('Axial distance [m]','fontsize',fnt);
            ylabel('Impact parameter [m]','fontsize',fnt);
            set(gca,'fontsize',fnt);
            set(gca,'xtick',[.077 .080 .083]);
            
            if sign_twin > 0
                saveas(fig1,[directory_save,'\den_int_contour_',num2str(shotnum_def),'_',num2str(shotnum_base),'POS.tif']);
                saveas(fig33,[directory_save,'\den_num_contour_',num2str(shotnum_def),'_',num2str(shotnum_base),'POS.tif']);
                saveas(fig401,[directory_save,'\ne_test_',num2str(shotnum_def),'_',num2str(shotnum_base),'POS.tif']);
            else
                saveas(fig1,[directory_save,'\den_int_contour_',num2str(shotnum_def),'_',num2str(shotnum_base),'.tif']);
                saveas(fig33,[directory_save,'\den_num_contour_',num2str(shotnum_def),'_',num2str(shotnum_base),'.tif']);
                saveas(fig401,[directory_save,'\ne_test_',num2str(shotnum_def),'_',num2str(shotnum_base),'.tif']);
            end
            
            
            
        else
        end
        
        %         check(shot_ind) = 1;
        %     catch
        %         check(shot_ind) = 0;
        %     end
        
    end

    
    save('interp_compare_nearest_mask','rad_l','rad_rec_l','N_prof_l','den_int_rec_l','den_num_rec_l','n_prof_l','ph_shift_l_short','ph_shift_l',...
        'rad_r','rad_rec_r','N_prof_r','den_int_rec_r','den_num_rec_r','n_prof_r','ph_shift_r_short','ph_shift_r',...
        'N_additional');
    
    
    for k = 1:length(cross_sect)
        den_int_split_l_tmp{k} = flipud(den_int_full(1:length(n_prof_l{k}),cross_sect(k)));
        den_int_split_l{k} = den_int_split_l_tmp{k}-den_int_split_l_tmp{k}(end)+N_edge_l(k);
        den_int_split_r_tmp{k} = den_int_full(length(n_prof_l{k})+1:end,cross_sect(k));
        den_int_split_r{k} = den_int_split_r_tmp{k}-den_int_split_r_tmp{k}(end)+N_edge_r(k);
    end
    % figure(444);
    % subplot(2,1,1);
    % plot(flipud(den_int_split_l{1}));
    % subplot(2,1,2);
    % plot(den_int_split_r{1});
    
    % figure(445);hold on;
    % subplot(2,1,1); hold on;
    % plot(N_after_l{1},'color','k');
    % plot(N_before_l{1},'color','m');
    % plot(den_int_rec_l{1},'color','g');
    % subplot(2,1,2); hold on;
    % plot(N_after_r{1},'color','k');
    % plot(N_before_r{1},'color','m');
    
    % k = 3;
    % figure(446);hold on;
    % suptitle(['k = ', num2str(k)]);
    % % subplot(2,2,1); hold on;
    % % plot(N_after_l{k},'color','k');
    % % plot(N_before_l{k},'color','m');
    % % plot(den_int_rec_l{k},'color','g');
    % % ylim([0 2e21]);
    % subplot(1,2,1); hold on;
    % plot(N_after_l{k},'color','k');
    % plot(interp1(rad_l{k},den_int_split_l{k},rad_rec_l{k}),'color','m','linewidth',2);
    % plot(interp1(rad_l{k},N_prof_l{k},rad_rec_l{k}),'color','g');
    % ylim([0 2e21]);
    %
    % % subplot(2,2,2); hold on;
    % % plot(N_after_r{k},'color','k');
    % % plot(N_before_r{k},'color','m');
    % % plot(den_int_rec_r{k},'color','g');
    % subplot(1,2,2); hold on;
    % plot(N_after_r{k},'color','k');
    % plot(interp1(rad_r{k},den_int_split_r{k},rad_rec_r{k}),'color','m');
    % plot(interp1(rad_r{k},N_prof_r{k},rad_rec_r{k}),'color','g');
    
    
    for k = 1:length(cross_sect)
        N_after_l{k} = den_int_rec_l{k}+N_additional(k);
        N_before_l{k} = interp1(rad_l{k},N_prof_l{k},rad_rec_l{k}); %using
        %     line-integrated density after Abel transform for comparison.
        %     N_before_l{k} = interp1(rad_l{k},den_int_split_l{k},rad_rec_l{k});
        
        N_error_l{k} = N_after_l{k} - N_before_l{k};
        if isnan(N_error_l{k}(end))
            N_error_l{k}(end)= N_error_l{k}(length(N_error_l{k})-1);
        else
        end
        
        %     if isnan(N_before_l{k}(end))
        %         N_before_l{k}(end) = N_before_l{k}(length(N_before_l{k})-1);
        %     else
        %     end
        %     N_error_l{k} = N_after_l{k} - N_before_l{k};
        
        n_error_out_l{k} = abel_invert_uncertainty(rad_rec_l{k},N_error_l{k});
        n_error_out_interp_l{k} = interp1(rad_rec_l{k},n_error_out_l{k},rad_l{k});
        
        %     n_error_high_l{k} = n_error_out_l{k};
        %     n_error_low_l{k} = n_error_out_l{k};
        %     ind_tmp_error{k} = find((den_num_rec_l{k}-n_error_low_l{k})<0);
        %     n_error_low_l{k}(ind_tmp_error{k}) = den_num_rec_l{k}(ind_tmp_error{k});
        
    end
    
    figure(1500);
    k = 1;
    subplot(3,1,k); hold on;
    title('Left side');
    h1500_2 = plot(rad_rec_l{k},N_before_l{k},'m','linewidth',3);
    h1500_1 = plot(rad_rec_l{k},N_after_l{k},'k','linewidth',3);
    % h1500_1 = errorbar(rad_rec_l{k},N_after_l{k},N_error_l{k},'color','k');
    xlim([min(rad_rec_l{k}) max(rad_rec_l{k})]);
    
    l1500 = legend([h1500_1 h1500_2],'N_after_l{k} = den_int_rec_l{k}+N_additional(k)','N_before_l{k} = interp1(rad_l{k},N_prof_l{k},rad_rec_l{k})');
    set(l1500,'Interpreter','none');
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},N_before_l{k},'m','linewidth',3);
        plot(rad_rec_l{k},N_after_l{k},'k','linewidth',3);
        % errorbar(rad_rec_l{k},N_after_l{k},N_error_l{k},'color','k');
        xlim([min(rad_rec_l{k}) max(rad_rec_l{k})]);
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},N_before_l{k},'m','linewidth',3);
        plot(rad_rec_l{k},N_after_l{k},'k','linewidth',3);
        %errorbar(rad_rec_l{k},N_after_l{k},N_error_l{k},'color','k');
        xlim([min(rad_rec_l{k}) max(rad_rec_l{k})]);
    else
    end
    
    for k = 1:length(cross_sect)
        N_after_r{k} = den_int_rec_r{k}+N_additional(k);
        N_before_r{k} = interp1(rad_r{k},N_prof_r{k},rad_rec_r{k}); %using
        %     line-integrated density after Abel transform for comparison.
        %     N_before_r{k} = interp1(rad_r{k},den_int_split_r{k},rad_rec_r{k});
        %     if isnan(N_before_r{k}(end))
        %         N_before_r{k}(end) = N_before_r{k}(length(N_before_r{k})-1);
        %     else
        %     end
        
        N_error_r{k} = N_after_r{k} - N_before_r{k};
        if isnan(N_error_r{k}(end))
            N_error_r{k}(end)= N_error_r{k}(length(N_error_r{k})-1);
        else
        end
        
        n_error_out_r{k} = abel_invert_uncertainty(rad_rec_r{k},N_error_r{k});
        n_error_out_interp_r{k} = interp1(rad_rec_r{k},n_error_out_r{k},rad_r{k});
        %     n_error_high_r{k} = n_error_out_r{k};
        %     n_error_low_r{k} = n_error_out_r{k};
        %     ind_tmp_error{k} = find((den_num_rec_r{k}-n_error_low_r{k})<0);
        %     n_error_low_r{k}(ind_tmp_error{k}) = den_num_rec_r{k}(ind_tmp_error{k})-1e22;
        
    end
    
    figure(1501);
    k = 1;
    subplot(3,1,k); hold on;
    title('Right side');
    h1501_2 = plot(rad_rec_r{k},N_before_r{k},'m','linewidth',3);
    h1501_1 = plot(rad_rec_r{k},N_after_r{k},'k','linewidth',3);
    % den_int_full_sect_r{k} = den_int_full(length(n_prof_l{k})+1:end,cross_sect(k));
    % plot(rad_r{k},den_int_full_sect_r{k}-min(den_int_full_sect_r{k}),'color','g');
    % h1501_1 = errorbar(rad_rec_r{k},N_after_r{k},N_error_r{k},'color','k');
    xlim([min(rad_rec_r{k}) max(rad_rec_r{k})]);
    l1501 = legend([h1501_1 h1501_2],'N_after_r{k} = den_int_rec_r{k}+N_additional(k)','N_before_r{k} = interp1(rad_r{k},N_prof_r{k},rad_rec_r{k})');
    set(l1501,'Interpreter','none');
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_rec_r{k},N_before_r{k},'m','linewidth',3);
        plot(rad_rec_r{k},N_after_r{k},'k','linewidth',3);
        %errorbar(rad_rec_r{k},N_after_r{k},N_error_r{k},'color','k');
        xlim([min(rad_rec_r{k}) max(rad_rec_r{k})]);
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_rec_r{k},N_before_r{k},'m','linewidth',3);
        plot(rad_rec_r{k},N_after_r{k},'k','linewidth',3);
        % errorbar(rad_rec_r{k},N_after_r{k},N_error_r{k},'color','k');
        xlim([min(rad_rec_r{k}) max(rad_rec_r{k})]);
    else
    end
    
    bar_spacing = 8;
    
    figure(1502);
    k = 1;
    subplot(3,1,k); hold on;
    plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
    shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
    % errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
    plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
    % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},[n_error_low_r{k} n_error_high_r{k}]','k');
    shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
    % errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
    
    ylim([0 4e23]);
    xlim([0 10e-3]);
    title(['z = 0.145 m'],'fontsize',fnt_2);
    ylabel('n_e [m^{-3}]','fontsize',fnt);
    set(gca,'fontsize',fnt_3);
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
        shadedErrorBar(rad_rec_l{k},den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
        plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
        shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5)
        % errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        title(['z = 0.150 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
        
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
        shadedErrorBar(rad_rec_l{k},den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
        plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
        shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5)
        % errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        title(['z = 0.155 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        xlabel('radius [m]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
    else
    end
    
    
    
    
    
    
    
    
    
    %plotting den_num_rec_s with error bars:
    figure(1503);
    k = 1;
    subplot(3,1,k); hold on;
    h1503_1 = plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
    % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
    bar_spacing = 5;
    errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
    h1503_2 = plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
    % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
    bar_spacing = 7;
    errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
    
    ylim([0 4e23]);
    xlim([0 12e-3]);
    title(['z = 0.145 m'],'fontsize',fnt_2);
    ylabel('n_e [m^{-3}]','fontsize',fnt);
    set(gca,'fontsize',fnt_3);
    l1503 = legend([h1503_1 h1503_2],'den_num_rec_l{k}','den_num_rec_r{k}');
    set(l1503,'Interpreter','none');
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k},den_num_rec_l{k},n_error_out_l{k},'r',.5)
        bar_spacing = 6;
        errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
        plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5)
        bar_spacing = 7;
        errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        title(['z = 0.150 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
        
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_rec_l{k},den_num_rec_l{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k},den_num_rec_l{k},n_error_out_l{k},'r',.5)
        bar_spacing = 3;
        errorbar(rad_rec_l{k}(1:bar_spacing:end),den_num_rec_l{k}(1:bar_spacing:end),n_error_out_l{k}(1:bar_spacing:end),'color','r');
        plot(rad_rec_r{k},den_num_rec_r{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5)
        bar_spacing = 7;
        errorbar(rad_rec_r{k}(1:bar_spacing:end),den_num_rec_r{k}(1:bar_spacing:end),n_error_out_r{k}(1:bar_spacing:end),'color','k');
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        title(['z = 0.155 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        xlabel('radius [m]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
    else
    end
    
    
    
    %plotting n_prof_s with error bars:
    figure(1504);
    k = 1;
    subplot(3,1,k); hold on;
    h1504_1 = plot(rad_l{k},n_prof_l{k},'r','linewidth',3);
    % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
    bar_spacing = 7;
    errorbar(rad_l{k}(1:bar_spacing:end),n_prof_l{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
    h1504_2 = plot(rad_r{k},n_prof_r{k},'k','linewidth',3);
    % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
    % bar_spacing = 10;
    errorbar(rad_r{k}(1:bar_spacing:end),n_prof_r{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
    
    ylim([0 4e23]);
    xlim([0 12e-3]);
    title(['z = 0.077 m'],'fontsize',fnt_2);
    ylabel('n_e [m^{-3}]','fontsize',fnt);
    set(gca,'fontsize',fnt_3);
    % l1504 = legend([h1504_1 h1504_2],'n_prof_l{k}','n_prof_r{k}');
    % set(l1504,'Interpreter','none');
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_l{k},n_prof_l{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % bar_spacing = 10;
        errorbar(rad_l{k}(1:bar_spacing:end),n_prof_l{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
        plot(rad_r{k},n_prof_r{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
        % bar_spacing = 10;
        errorbar(rad_r{k}(1:bar_spacing:end),n_prof_r{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
        
        ylim([0 4e23]);
        xlim([0 12e-3]);
        title(['z = 0.080 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
        
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_l{k},n_prof_l{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % bar_spacing = 10;
        errorbar(rad_l{k}(1:bar_spacing:end),n_prof_l{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
        plot(rad_r{k},n_prof_r{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
        % bar_spacing = 10;
        errorbar(rad_r{k}(1:bar_spacing:end),n_prof_r{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
        
        
        ylim([0 4e23]);
        xlim([0 12e-3]);
        title(['z = 0.083 m'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        xlabel('radius [m]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
    else
    end
    
    
    
    % computing linear density for comparison with ZaP:
    k = 1;
    rad_int_ind = find(rad_l{k}>.0025,1,'first');
    N1_tmp = cumtrapz(rad_l{k}(1:rad_int_ind),2*pi*den_num_axial_l_plt{k}(1:rad_int_ind).*rad_l{k}(1:rad_int_ind)');
    N1_end = N1_tmp(end)
    
    % plotting the density profiles with axial variation included:
    figure(1505);
    k = 1;
    subplot(3,1,k); hold on;
    h1504_1 = plot(rad_l{k},den_num_axial_l_plt{k},'r','linewidth',3);
    % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
    bar_spacing = 7;
    errorbar(rad_l{k}(1:bar_spacing:end),den_num_axial_l_plt{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
    h1504_2 = plot(rad_r{k},den_num_axial_r_plt{k},'k','linewidth',3);
    % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
    % bar_spacing = 10;
    errorbar(rad_r{k}(1:bar_spacing:end),den_num_axial_r_plt{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
    
    ylim([0 4e23]);
    xlim([0 10e-3]);
    % title(['z = 0.077 m'],'fontsize',fnt_2);
    % title(['(a)'],'fontsize',fnt_2);
    ylabel('n_e [m^{-3}]','fontsize',fnt);
    % xlabel([{'radius [m]'},{'(a)'}],'fontsize',fnt);
    set(gca,'fontsize',fnt_3);
    % l1504 = legend([h1504_1 h1504_2],'n_prof_l{k}','n_prof_r{k}');
    % set(l1504,'Interpreter','none');
    
    if length(cross_sect) ==3
        k = 2;
        subplot(3,1,k); hold on;
        plot(rad_l{k},den_num_axial_l_plt{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % bar_spacing = 10;
        errorbar(rad_l{k}(1:bar_spacing:end),den_num_axial_l_plt{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
        plot(rad_r{k},den_num_axial_r_plt{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
        % bar_spacing = 10;
        errorbar(rad_r{k}(1:bar_spacing:end),den_num_axial_r_plt{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        % title(['z = 0.080 m'],'fontsize',fnt_2);
        % title(['(b)'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        % xlabel([{'radius [m]'},{'(b)'}],'fontsize',fnt);
        set(gca,'fontsize',fnt_3);
        
        k = 3;
        subplot(3,1,k); hold on;
        plot(rad_l{k},den_num_axial_l_plt{k},'r','linewidth',3);
        % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
        % bar_spacing = 10;
        errorbar(rad_l{k}(1:bar_spacing:end),den_num_axial_l_plt{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',3);
        plot(rad_r{k},den_num_axial_r_plt{k},'k','linewidth',3);
        % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
        % bar_spacing = 10;
        errorbar(rad_r{k}(1:bar_spacing:end),den_num_axial_r_plt{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',3);
        
        
        ylim([0 4e23]);
        xlim([0 10e-3]);
        % title(['z = 0.083 m'],'fontsize',fnt_2);
        % title(['(c)'],'fontsize',fnt_2);
        ylabel('n_e [m^{-3}]','fontsize',fnt);
        xlabel('radius [m]','fontsize',fnt);
        set(gca,'fontsize',fnt_3);
    else
    end
    
    figure(1506); hold on;
    fnt_3 = 24;
    k = 1;
    % subplot(3,1,k); hold on;
    wid = 7;
    h1504_1 = plot(rad_l{k},den_num_axial_l_plt{k},'r','linewidth',wid);
    % shadedErrorBar(rad_rec_l{k}',den_num_rec_l{k},n_error_out_l{k},'r',.5)
    bar_spacing = 7;
    errorbar(rad_l{k}(1:bar_spacing:end),den_num_axial_l_plt{k}(1:bar_spacing:end),n_error_out_interp_l{k}(1:bar_spacing:end),'color','r','linestyle','none','linewidth',wid);
    h1504_2 = plot(rad_r{k},den_num_axial_r_plt{k},'k','linewidth',wid);
    % shadedErrorBar(rad_rec_r{k},den_num_rec_r{k},n_error_out_r{k},'k',.5);
    % bar_spacing = 10;
    errorbar(rad_r{k}(1:bar_spacing:end),den_num_axial_r_plt{k}(1:bar_spacing:end),n_error_out_interp_r{k}(1:bar_spacing:end),'color','k','linestyle','none','linewidth',wid);
    
    ylim([0 4e23]);
    xlim([0 10e-3]);
    % title(['z = 0.077 m'],'fontsize',fnt_2);
    % title(['(a)'],'fontsize',fnt_2);
    ylabel('n_e [m^{-3}]','fontsize',fnt);
    % xlabel([{'radius [m]'},{'(a)'}],'fontsize',fnt);
    set(gca,'fontsize',fnt_3);
    % l1504 = legend([h1504_1 h1504_2],'n_prof_l{k}','n_prof_r{k}');
    % set(l1504,'Interpreter','none');
    xlabel('radius [m]','fontsize',fnt_3);
    set(gca,'xtick',[0 0.0025 0.005 0.0075 0.01]);
    
    figure(112);hold on;
    suptitle('phase');
    subplot(1,2,1);hold on;
    plot(unwrapped(1,:));
    title('left');
    subplot(1,2,2);hold on;
    plot(unwrapped(end,:));
    title('right');
    
    figure(113);hold on;
    suptitle('Ne');
    subplot(1,2,1);hold on;
    plot(den_int_full(1,:));
    title('left');
    subplot(1,2,2);hold on;
    plot(den_int_full(end,:));
    title('right');
    
    figure(114);hold on;
    suptitle('Ne');
    subplot(1,2,1);hold on;
    plot(den_num_full_axial(1,:));
    title('left');
    subplot(1,2,2);hold on;
    plot(den_num_full_axial(end,:));
    title('right');
    
    figure(115);
    waterfall(x_twin_red_adj,y_twin_red_adj,den_num_full_axial);
    figure(116);
    waterfall(y_twin_red_adj,x_twin_red_adj,den_num_full_axial');
    
    
    
    
    % plot the abel inverted profiles that Ray wanted to see:
    
    
    
    for k = 1:length(cross_sect)
        den_int_ray_l{k} = flipud(den_int_full(1:length(den_num_l{cross_sect(k)}),cross_sect(k)));
        den_int_ray_r{k} = den_int_full((length(den_num_l{cross_sect(k)})+1):end,cross_sect(k));
    end
    
    figure(991);hold on;
    k = 1;
    plot(den_int_ray_l{k});
    plot(den_int_ray_r{k});
    
    for k = 1:length(cross_sect)
        den_num_ray_l{k} = abel_invert(rad_l{k},den_int_ray_l{k}-min(den_int_ray_l{k}));
        den_num_ray_r{k} = abel_invert(rad_r{k},den_int_ray_r{k}-min(den_int_ray_r{k}));
        den_num_ray_l_2{k} = abel_invert(rad_l{k},den_int_ray_l{k});
        den_num_ray_r_2{k} = abel_invert(rad_r{k},den_int_ray_r{k});
    end
    
    figure(992); hold on;
    k = 1;
    plot(rad_l{k},den_num_ray_l{k},'r:');
    plot(rad_r{k},den_num_ray_r{k},'k:');
    plot(rad_l{k},den_num_ray_l_2{k},'r--');
    plot(rad_r{k},den_num_ray_r_2{k},'k--');
    plot(rad_l{k},den_num_axial_l_plt{k},'r','linewidth',3);
    plot(rad_r{k},den_num_axial_r_plt{k},'k','linewidth',3);
    
    for k = 1:length(cross_sect)
        den_int_ray_abel_l{k} = abel(rad_l{k},den_num_axial_l_plt{k}');
        den_int_ray_abel_r{k} = abel(rad_r{k},den_num_axial_r_plt{k}');
        den_int_ray_abel_l_2{k} = abel(rad_l{k},den_num_axial_l_plt{k}'-min(den_num_axial_l_plt{k}'));
        den_int_ray_abel_r_2{k} = abel(rad_r{k},den_num_axial_r_plt{k}'-min(den_num_axial_r_plt{k}'));
    end
    
    figure(993); hold on;
    k = 1;
    plot(rad_l{k},den_int_ray_l{k},'r','linewidth',3);
    plot(rad_r{k},den_int_ray_r{k},'k','linewidth',3);
    plot(rad_l{k},den_int_ray_abel_l{k},'r:');
    plot(rad_r{k},den_int_ray_abel_r{k},'k:');
    plot(rad_l{k},den_int_ray_abel_l_2{k}+7.8114e21,'g','linewidth',3);
    plot(rad_r{k},den_int_ray_abel_r_2{k}+1.5939e22,'g','linewidth',3);
    
    test_prof_111_l = den_int_ray_abel_l_2{k}+7.8114e21;
    test_prof_111_r = den_int_ray_abel_r_2{k}+1.5939e22;
    %
    % test_prof_111_trans_l = abel(rad_l{k},test_prof_111_l'-min(test_prof_111_l));
    % test_prof_111_trans_r = abel(rad_r{k},test_prof_111_r'-min(test_prof_111_r));
    
    test_prof_111_trans_l = abel_invert(rad_l{k},test_prof_111_l-min(test_prof_111_l))+1.924e22;
    test_prof_111_trans_r = abel_invert(rad_r{k},test_prof_111_r-min(test_prof_111_r))+3.924e22;
    
    test_prof_111_abel_l = abel(rad_l{k},test_prof_111_trans_l'-min(test_prof_111_trans_l));
    test_prof_111_abel_r = abel(rad_r{k},test_prof_111_trans_r'-min(test_prof_111_trans_r));
    
    
    figure(123423);
    subplot(3,1,1);hold on;
    plot(rad_l{k},test_prof_111_l);
    plot(rad_r{k},test_prof_111_r);
    plot(rad_l{k},den_int_ray_l{k},'color','k');
    plot(rad_r{k},den_int_ray_r{k},'color','k');
    
    subplot(3,1,2);hold on;
    plot(rad_l{k},test_prof_111_trans_l);
    plot(rad_r{k},test_prof_111_trans_r);
    subplot(3,1,3);hold on;
    plot(rad_l{k},test_prof_111_abel_l);
    plot(rad_r{k},test_prof_111_abel_r);
    
    k = 1;
    for i = 1:5
        den_int_test_plt{i} = abel(rad_l{k},den_num_axial_l_plt{k}'-min(den_num_axial_l_plt{k}')+i*5e22-5e22);    
    end
    
    figure(995);hold on;
    plot(rad_l{k},den_int_test_plt{1},'linewidth',3);
    plot(rad_l{k},den_int_test_plt{2});
    plot(rad_l{k},den_int_test_plt{3});
    plot(rad_l{k},den_int_test_plt{4});
    plot(rad_l{k},den_int_test_plt{5});