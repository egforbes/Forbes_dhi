% DHI_recon_invert_all.m 
% Author:  Michael P. Ross
% Date = Sept. 6, 2016
% Description:  This code takes the shot numbers and reconstruction 
% settings specified in the shotnumber_inputs function and performs a 
% Fresnel transform reconstruction on the data.  
%
% Then, it proceeds to Abel invert the line-integrated density data to
% obtain number density.  The code follows the procedure described in 
% Sec. 3.6.2 to minimize the difference in the resulting sides of the
% number density profile.
%
% The code saves some data as a .mat file for use in other analysis codes.
% It saves images of the reconstructed phase, line-integrated density, and
% number density to the holography network drive.  Finally, it uses the
% function MDS_DHI_put.m to store the density data alongside reconstruction
% settings in ZaP-HD's MDSplus data tree.

clear all; clc;
close all;

%% User inputs, file directories, and constants:

% Importing the desired shot numbers:
[shots,d,xmin,xmax,ymin,ymax,sign_twin] = shotnumber_inputs(2,1);
% Fraction along z axis for each cross-section:
%%
cross_sect_frac = [.25 0.5 0.75]; 
% Axial measurement location:
z_loc = 0.10; % [m]

% Creating a directory for .mat file and image storage:
date_ref = num2str(shots(1,1));
date = date_ref(1:6);
directory_save = ['M:\Users\Eleanor_Forbes\Reconstructions\',date,'\inverted'];
mkdir(['M:\Users\Eleanor_Forbes\Reconstructions\',date],'inverted');

% Declare constants and ZaP-HD plasma parameters:
e_0 = 8.85e-12; % Permitivity of free space
mu_0 = 1.257e-6; % Permeability of free space
k = 1.38e-23; % Boltzmann's constant in SI units
q = 1.6e-19; % electron charge [C]
c = 3e8; % speed of light [m/s]
mi = 1.67e-27; % proton mass [kg]
me = 9.12e-31; % electron mass [kg]
lambda = 532e-9; % laser wavelength [m]
R_electrode = 8*2.54/100;

%% 2.) Fresnel transform reconstruction and Abel inversion:
% Iterate through desired shots:
for shot_ind = 1:size(shots,1)
    % Close figures and clear data on each iteration for memory management:
    close all;
    clearvars -except h theta d_fringe step_sm Bw V ni ne Ti Te T r ...
        R_electrode lambda me mi c q k mu_0 e_0 z_loc cross_sect_frac ...
        shots d xmin xmax ymin ymax sign_twin shot_ind directory_save date
    
    shotnum_base = shots(shot_ind,2);
    shotnum_def = shots(shot_ind,1);
    
    %% Fresnel transform reconsruction:
    [unwrapped,x_twin,y_twin,phase_diff,phase_x,phase_y] = ...
        holographic_reconstruction(shotnum_base,shotnum_def,xmin,...
        xmax,ymin,ymax,d);
    % Plot and save the reconstructed phase:
    fig999 = figure(999);
    imagesc(phase_x,phase_y,phase_diff);
%    saveas(fig999,[directory_save,'\phase_map_',num2str(...
%         shotnum_def),'_',num2str(shotnum_base),'POS.tif']);
    
    % Define the indices for the selected cross sections:
    cross_sect = round(cross_sect_frac*size(unwrapped,2));
    
    % Converting the phase to the 2-D line-integrated density:
    % Note that the sign of den_int_full here depends on the optical
    % set up (e.g. the twin images can flip if the scene and reference
    % beams flip orientation).
    den_int_full = sign_twin*[unwrapped/((-q^2/(4*pi*c^2*me*e_0))...
        *lambda)]';
    den_int_full = den_int_full - min(min(den_int_full));
    
    % Reducing the size of x_twin and y_twin to match unwrapped
    % (unwrapped is reduced in size because of boxcar smoothing):
    del_y = y_twin(2) - y_twin(1);
    del_x = x_twin(2) - x_twin(1);
    y_twin_red = [0:del_y:del_y*(size(den_int_full,2)-1)];
    x_twin_red = [0:del_x:del_x*(size(den_int_full,1)-1)];
    
    % Adjusting the zero value of x_twin and y_twin so the origin is 
    % centered.  Can adjust here for when the hologram is not centered at 
    % 0 impact param.
    y_twin_red_adj = (y_twin_red - y_twin_red(end)/2);
    x_twin_red_adj = x_twin_red - x_twin_red(end)/2+z_loc;
  
    
    % Save .mat file with line-integrated data for further analysis:
%     save(['M:\Users\Eleanor_Forbes\Reconstructions\',date,'\Ne_',num2str(shotnum_def),...
%         '.mat'],'x_twin_red_adj','y_twin_red_adj','den_int_full');
    
    % Identify number of cross sections:
    num_of_cross_sect = size(unwrapped,2);
    
    % Identify number of centroids:
    num_of_centroids = 25;
    
    % Identify number of n_edge values:
    num_of_n_edge = 25;
    n_edge_delta = 1e20;
    
    %% Peform the Abel inversion:
    
    % Iterate through cross sections:
%     for k = 1:num_of_cross_sect

for k = 1:num_of_cross_sect
        
        % Identify the max Ne value and its index:
        den_int = den_int_full(k,:);
        [Ne_max(k),ind_Ne_max(k)] = max(den_int);
        
        % Initialize centroid_ind_tmp to enable determination of 
        % centroid_abs (the centroid locations in absolute pixel-space):
        centroid_ind_tmp = ind_Ne_max(k)-num_of_centroids;
        
        % Reset the centroid_iterations if it was reduced to accomodate
        % centroids close to the bounds of the data:
        centroids_iterations = num_of_centroids;
       
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
            
       
            %%
            % Iterate through centroids:
            for i = 1:centroids_iterations
                
                % Selecting the left side of the Ne profile and flipping it
                Ne_1d_left_temp = fliplr(den_int(1:centroid_ind_tmp));
                
                % Selecting the right side of the Ne profile:
                Ne_1d_right_temp = den_int(centroid_ind_tmp+1:end);
%                 figure(1)
%                 plot(Ne_1d_left{k,i})
%                 hold on
%                 plot(Ne_1d_right{k,i})
%                 disp('Press a key')
%                 pause      

 % ELLIE CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ne_min_l,ind_Ne_min_l] = min(Ne_1d_left_temp);
[Ne_min_r,ind_Ne_min_r] = min(Ne_1d_right_temp);

if ind_Ne_min_l ~= length(Ne_1d_left_temp)
    truncate_l{k,i} = length(Ne_1d_left_temp) - ind_Ne_min_l;
    Ne_1d_left{k,i} = Ne_1d_left_temp(1:ind_Ne_min_l);
else
    truncate_l{k,i} = 0;
    Ne_1d_left{k,i} = Ne_1d_left_temp;
end
if ind_Ne_min_r ~= length(Ne_1d_right_temp)
    truncate_r{k,i} = length(Ne_1d_right_temp) - ind_Ne_min_r;
    Ne_1d_right{k,i} = Ne_1d_right_temp(1:ind_Ne_min_r);
else
    truncate_r{k,i} = 0;
    Ne_1d_right{k,i} = Ne_1d_right_temp;
end

% END ELLIE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                
                % Iterate the index to prepare the code to move to the next
                % centroid value on the next loop iteration:
                centroid_ind_track(k,i) = centroid_ind_tmp;
                centroid_ind_tmp = centroid_ind_tmp + 1;
                
                % Conductig the Abel inversion to compute the ne profile on
                % each side of the assumed centroid:
                ne_1d_left{k,i} = abel_invert(y_twin,[Ne_1d_left{k,i}-...
                    min(Ne_1d_left{k,i})]'); 
                ne_1d_right{k,i} = abel_invert(y_twin,[Ne_1d_right{k,i}-...
                    min(Ne_1d_right{k,i})]');
%                 
%                 figure(2)
%                 plot(ne_1d_left{k,i})
%                 hold on
%                 plot(ne_1d_right{k,i})
%                 disp('Press a key')
%                 pause
                
                
                % When the centroid does not divide the Ne profile in half,
                % the longer of the two sides must be truncated to match 
                % the length of the shorter side.  The following case 
                % statements conduct this truncation:
                len_left(k,i) = length(Ne_1d_left{k,i});
                len_right(k,i) = length(Ne_1d_right{k,i});
                
                if len_left(k,i)>len_right(k,i)
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i}(1:...
                        len_right(k,i));
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i};
                elseif len_left(k,i)<len_right(k,i)
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i};
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i}(1:...
                        len_left(k,i));
                else
                    ne_1d_left_trunc{k,i} = ne_1d_left{k,i};
                    ne_1d_right_trunc{k,i} = ne_1d_right{k,i};
                end
                
                % Identify the nominal n_edge_0 value as the density in the
                % longer profile at the radius of the edge of the short 
                % profile:
                if len_left(k,i)>len_right(k,i) && len_right(k,i)~=0
                    n_edge_0(k,i) = ne_1d_left_trunc{k,i}(len_right(k,i));
                elseif len_left(k,i)<len_right(k,i) && len_left(k,i)~=0
                    n_edge_0(k,i) = ne_1d_right_trunc{k,i}(len_left(k,i));
                else
                    n_edge_0(k,i) = 0;
                end
                
                % Initialize n_edge_tmp to enable determination of 
                % centroid_abs (the centroid locations in absolute 
                % pixel-space):
                n_edge_ind_tmp = -num_of_n_edge;
                
                % Iterate through n_edge values:
                for j = 1:2*num_of_n_edge
                    % populating list of all n_edge attempted:
                    n_edge_val(j) = n_edge_ind_tmp*n_edge_delta;
                    
                    % Iterating through n_edge values for the shorter 
                    % profile:
                    if len_left(k,i)>len_right(k,i)
                        ne_1d_left_compare{k,i,j} = ne_1d_left_trunc{k,i};
                        ne_1d_right_compare{k,i,j} = ...
                            ne_1d_right_trunc{k,i}+...
                            n_edge_delta*n_edge_ind_tmp;
                    elseif len_left(k,i)<len_right(k,i)
                        ne_1d_left_compare{k,i,j} = ...
                            ne_1d_left_trunc{k,i}+...
                            n_edge_delta*n_edge_ind_tmp;
                        ne_1d_right_compare{k,i,j} = ...
                            ne_1d_right_trunc{k,i};
                    else
                        if ne_1d_left_trunc{k,i}(end) > ...
                                ne_1d_right_trunc{k,i}(end)
                            ne_1d_left_compare{k,i,j} = ...
                                ne_1d_left_trunc{k,i}+...
                                n_edge_delta*n_edge_ind_tmp;
                            ne_1d_right_compare{k,i,j} = ...
                                ne_1d_right_trunc{k,i};
                        else
                            ne_1d_left_compare{k,i,j} = ...
                                ne_1d_left_trunc{k,i};
                            ne_1d_right_compare{k,i,j} = ...
                                ne_1d_right_trunc{k,i}+...
                                n_edge_delta*n_edge_ind_tmp;
                        end
                    end
                    
                    % Iterate the index to prepare the code to move to the 
                    % next n_edge value on the next loop iteration:
                    n_edge_ind_track(k,i,j) = n_edge_ind_tmp;
                    n_edge_ind_tmp = n_edge_ind_tmp + 1;
                    
                    % Computing the L2 norm of the difference in the sides 
                    % of the ne profile:
                    dens_diff{k,i,j} = abs(ne_1d_left_compare{k,i,j}(1:2)-...
                        ne_1d_right_compare{k,i,j}(1:2));
                    norm_dens(k,i,j) = norm(dens_diff{k,i,j},2)/...
                        (length(dens_diff{k,i,j}(1:2)));
%                     dens_diff{k,i,j} = abs(ne_1d_left_compare{k,i,j}-...
%                         ne_1d_right_compare{k,i,j});
%                     norm_dens(k,i,j) = norm(dens_diff{k,i,j},2)/...
%                         (length(dens_diff{k,i,j})); 
%                     figure(4)
%                         plot(ne_1d_left_compare{k,i,j})
%                         hold on
%                         plot(ne_1d_right_compare{k,i,j})
                end
                
                % Changing all zero-valued elements in norm_dens to NaNs 
                % (allows the code to ignore them when searching for the 
                % minimum)
                ind_zeros = find(norm_dens(k,i,:)==0);
                norm_dens(k,i,ind_zeros) = NaN;
                
                % Find the minimum L2 norm value for each combination of 
                % cross-section and centroid.  (n_edge_intermediate records
                % the index location of n_edge yielding minimum L2 norm)  
                [norm_intermediate(k,i),n_edge_intermediate(k,i)] = ...
                    min(norm_dens(k,i,:),[],3); 
            end
            
            % Find the minimum L2 norm value for each cross-section.  
            % (centroid records the relative index location of the centroid
            % location yielding minimum L2 norm)
            [norm_min_val(k),centroid(k)] = min(norm_intermediate(k,:)...
                ,[],2); 
            n_edge(k) = n_edge_intermediate(k,centroid(k));
%             pause
%             disp(centroid(k))
%             disp('Press a key')
            
            % Assigning the final left and right number density profiles:
            if centroid(k) == 0 % do not adjust for axial density variation
                %for non-inverted cross-sections
                
            elseif len_left(k,centroid(k))>len_right(k,centroid(k))
                den_num_l{k} = ne_1d_left{k,centroid(k)};
                den_num_r{k} = ne_1d_right{k,centroid(k)}+n_edge_delta*...
                    (n_edge(k)-num_of_n_edge);
                den_int_edge_axial(k) = den_int_full(k,1);
                
            elseif len_left(k,centroid(k))<len_right(k,centroid(k))
                den_num_l{k} = ne_1d_left{k,centroid(k)}+n_edge_delta*...
                    (n_edge(k)-num_of_n_edge);
                den_num_r{k} = ne_1d_right{k,centroid(k)};
                den_int_edge_axial(k) = den_int_full(k,end);
                
            else
                if ne_1d_left_trunc{k,i}(end) > ne_1d_right_trunc{k,i}(end)
                    den_num_l{k} = ne_1d_left_trunc{k,centroid(k)}+...
                        n_edge_delta*(n_edge(k)-num_of_n_edge);
                    den_num_r{k} = ne_1d_right_trunc{k,centroid(k)};
                    den_int_edge_axial(k) = den_int_full(k,end);
                    
                else
                    den_num_l{k} = ne_1d_left_trunc{k,centroid(k)};
                    den_num_r{k} = ne_1d_right_trunc{k,centroid(k)}+...
                        n_edge_delta*(n_edge(k)-num_of_n_edge);
                    den_int_edge_axial(k) = den_int_full(k,end);
                    
                end
                
                %Ellie code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if truncate_l{k,centroid(k)} ~= 0
                    den_num_l{k} = [den_num_l{k}; NaN*ones(truncate_l{k,centroid(k)},1)];
                end
                if truncate_r{k,centroid(k)} ~= 0
                    den_num_r{k} = [den_num_r{k}; NaN*ones(truncate_r{k,centroid(k)},1)];
                end
                %End Ellie Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            % Computing the centroid indices in absolute pixel coordinates:
            if centroid_ind_track(k,centroid(k)) == 0
                centroid_abs(k) = size(unwrapped,1);
            else
                centroid_abs(k) = centroid_ind_track(k,centroid(k));
            end
        end
        

                
        
    end
    
    if quit2nextShot == 0 % for debugging
        
        
        
 %% Compiling the inverted number density data to form the 
        % den_num_full matrix:
        den_num_full = zeros(size(den_int_full));
        for k = 1:num_of_cross_sect
            if centroid_abs(k) == 1 || centroid_abs(k) == size(unwrapped,1)
                den_num_full(k,:) = zeros(size(unwrapped,1),1);
            else
                den_num_full(k,(1+centroid_abs(k)-length(den_num_l{k})):...
                    centroid_abs(k)) = flipud(den_num_l{k});
                den_num_full(k,(1+centroid_abs(k)):(centroid_abs(k)+...
                    length(den_num_r{k}))) = den_num_r{k};
            end
        end
        
        %% Plotting the data:
        x_0 = x_twin_red_adj(1);
        y_0 = y_twin_red_adj(1);
        
        den_num_full = den_num_full';
        den_int_full = den_int_full';
        
        fnt = 24;
        
        % Plotting Ne (den_int_full):
        fig1 = figure(1);hold on;
        h11 = pcolor(x_twin_red_adj,y_twin_red_adj,den_int_full);
        set(h11,'edgecolor','none');
        colormap jet
        set(fig1,'position',[1290 520 560 420]);
        colorbar;
        xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
        ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
        title('Line-integrated electron density, N_e [m^{-2}]',...
            'fontsize',fnt);
        xlabel('Axial distance [m]','fontsize',fnt);
        ylabel('Impact parameter [m]','fontsize',fnt);
        set(gca,'fontsize',fnt);
        set(gca,'xtick',[0.145 0.15 0.155]);
        
        % Plotting number density not adjusted for axial variation:
        fig3 = figure(3);hold on;
        h31 = pcolor(x_twin_red_adj,y_twin_red_adj,den_num_full);
        set(h31,'edgecolor','none');
        h32 = scatter(x_twin_red_adj,y_twin_red_adj(centroid_abs),'k','.');
        set(h32,'sizedata',1000);
        colormap jet
        set(fig3,'position',[1867         519         560         420]);
        colorbar;
        xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
        ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
        line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],...
            [y_twin_red_adj(1) y_twin_red_adj(end)],'color','k');
        line([cross_sect(1)*del_x+x_0 cross_sect(1)*del_x+x_0],...
            [y_twin_red_adj(1) y_twin_red_adj(end)],'color','k',...
            'linewidth',3);
        if length(cross_sect)==3
            line([cross_sect(2)*del_x+x_0 cross_sect(2)*del_x+x_0],...
                [y_twin_red_adj(1) y_twin_red_adj(end)],'color','k',...
                'linewidth',3);%,'linestyle',':');
            line([cross_sect(3)*del_x+x_0 cross_sect(3)*del_x+x_0],...
                [y_twin_red_adj(1) y_twin_red_adj(end)],'color','k',...
                'linewidth',3);%,'linestyle','--');
        else
        end
        title('Electron number density, n_e [m^{-3}]','fontsize',fnt);
        xlabel('Axial distance [m]','fontsize',fnt);
        ylabel('Impact parameter [m]','fontsize',fnt);
        set(gca,'fontsize',fnt);
        set(gca,'xtick',[0.145 0.15 0.155]);
        
        %% Adjusting number density to account for axial variation of phase
        if exist('centroid') % do not adjust for axial density variation 
            % for non-inverted cross-sections
            
            for k = 1:num_of_cross_sect
                
                if centroid_abs(k) == 1
                    y_edge_l{k} = 0:del_y:del_y;
                    y_edge_r{k} = 0:del_y:size(unwrapped,1)-1;
                    y_count(k) = max(y_edge_l{k});
                    
                    x_edge_l(k) = 2*sqrt(R_electrode^2-max(y_edge_l{k})^2);
                    x_edge_r(k) = 2*sqrt(R_electrode^2-max(y_edge_r{k})^2);
                    
                elseif centroid_abs(k) == size(unwrapped,1)
                    y_edge_l{k} = 0:del_y:(size(unwrapped,1)-1)*del_y;
                    y_edge_r{k} = 0:del_y:del_y;
                    y_count(k) = max(y_edge_l{k});
                    
                    x_edge_l(k) = 2*sqrt(R_electrode^2-max(y_edge_l{k})^2);
                    x_edge_r(k) = 2*sqrt(R_electrode^2-max(y_edge_r{k})^2);
                    
                else
                    y_edge_l{k} = 0:del_y:(length(den_num_l{k})-1)*del_y;  
                    y_edge_r{k} = 0:del_y:(length(den_num_r{k})-1)*del_y;
                    
                    y_count(k) = max(y_edge_l{k});
                    x_edge_l(k) = 2*sqrt(R_electrode^2-max(y_edge_l{k})^2);
                    x_edge_r(k) = 2*sqrt(R_electrode^2-max(y_edge_r{k})^2);
                end
                
                if centroid(k) == 0 % do not adjust for axial density 
                    %variation for non-inverted cross-sections
                    
                elseif len_left(k,centroid(k))>len_right(k,centroid(k))%%
                    den_num_edge_axial_l(k) = den_int_edge_axial(k)/...
                        x_edge_l(k);
                    den_num_edge_axial_r(k) = 0;
                    
                    den_num_axial_l{k} = den_num_l{k} + ...
                        den_num_edge_axial_l(k);
                    den_num_axial_r{k} = den_num_r{k} + ...
                        den_num_edge_axial_l(k);
                    
                elseif len_left(k,centroid(k))<len_right(k,centroid(k))
                    den_num_edge_axial_l(k) = 0;
                    den_num_edge_axial_r(k) = den_int_edge_axial(k)/...
                        x_edge_l(k);
                    
                    den_num_axial_l{k} = den_num_l{k} + ...
                        den_num_edge_axial_r(k);
                    den_num_axial_r{k} = den_num_r{k} + ...
                        den_num_edge_axial_r(k);
                    
                else
                    den_num_edge_axial_l(k) = 0;
                    den_num_edge_axial_r(k) = 0;
                    
                    den_num_axial_l{k} = den_num_l{k};
                    den_num_axial_r{k} = den_num_r{k};
                end
            end
            
            % Compiling the adjusted number density into a single matrix:
            for k = 1:num_of_cross_sect
                
                if centroid(k) ==0
                    den_num_full_axial(k,1:size(unwrapped,1)) = ...
                        zeros(size(unwrapped,1),1);
                else
                    den_num_full_axial(k,(1+centroid_abs(k)-...
                        length(den_num_axial_l{k})):centroid_abs(k)) = ...
                        flipud(den_num_axial_l{k});
                    den_num_full_axial(k,(1+centroid_abs(k)):...
                        (centroid_abs(k)+length(den_num_axial_r{k}))) = ...
                        den_num_axial_r{k};
                end
            end
            den_num_full_axial = den_num_full_axial';
            
            %% Plotting the number density adjusted for axial variation:
            fig33 = figure(33);hold on; fnt = 24;
            h31 = pcolor(x_twin_red_adj,y_twin_red_adj,den_num_full_axial);
            set(h31,'edgecolor','none');
            h32 = scatter(x_twin_red_adj,y_twin_red_adj(centroid_abs),...
                'k','.');
            set(h32,'sizedata',1000);
            colormap jet
            set(fig33,'position',[1867 519 560 420]);
            colorbar;
            xlim([x_twin_red_adj(1) x_twin_red_adj(end)]);
            ylim([y_twin_red_adj(1) y_twin_red_adj(end)]);
            title('Electron number density, n_e [m^{-3}]','fontsize',fnt);
            xlabel('Axial distance [m]','fontsize',fnt);
            ylabel('Impact parameter [m]','fontsize',fnt);
            set(gca,'fontsize',fnt);
            set(gca,'xtick',[0.145 0.15 0.155]);
            legend(h32,'centroid');
            
            % Saving figures as images:
            if sign_twin > 0
                saveas(fig1,[directory_save,'\den_int_contour_',...
                    num2str(shotnum_def),'_',num2str(shotnum_base),...
                    'POS.tif']);
                saveas(fig33,[directory_save,'\den_num_contour_',...
                    num2str(shotnum_def),'_',num2str(shotnum_base),...
                    'POS.tif']);
            else
                saveas(fig1,[directory_save,'\den_int_contour_',...
                    num2str(shotnum_def),'_',num2str(shotnum_base),...
                    '.tif']);
                saveas(fig33,[directory_save,'\den_num_contour_',...
                    num2str(shotnum_def),'_',num2str(shotnum_base),...
                    '.tif']);
            end
            
            for k = 1:num_of_cross_sect
                rad_l_tmp{k} = 0:del_y:(length(den_num_l{k})-1)*del_y;  
                % del_y is in empirical reconstruction coords.
                rad_r_tmp{k} = 0:del_y:(length(den_num_r{k})-1)*del_y; 
            end
            for k = 1:length(cross_sect)
                rad_l{k} = rad_l_tmp{cross_sect(k)};
                rad_r{k} = rad_r_tmp{cross_sect(k)};
            end
            
        else
        end
                  
    else
    end
    
    %% Writing data to the ZaP-HD MDSplus data tree:
    %commented until MDS error can be resolved EF
%     MDS_DHI_put(shotnum_def,shotnum_base,den_int_full,den_num_full,...
%         centroid_abs,d,xmin,xmax,ymin,ymax,sign_twin,...
%         x_twin_red_adj,y_twin_red_adj)
end     

