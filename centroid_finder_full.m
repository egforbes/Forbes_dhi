%New centroid finding routine for entire 2D profile (small)
clear all; close all; clc;

%load the data
load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\backup_good','\Ne_',num2str(171005023),'_backup.mat'))

for jcol = 1:size(den_int_full,1)
    den_int = den_int_full(jcol,:); %pick out a line through the pinch
    [Ne_max, ind_Ne_max] = max(den_int); %find peak density location
    num_of_centroids = 6; %number of centroid guesses for each loop
    
    n_tries = 16; %number of times to add/subtract density from one profile. MUST BE EVEN
    
    centroid_ind_tmp = ind_Ne_max-num_of_centroids; %Starting location of the centroid
    
    for jcen = 1:2*num_of_centroids
        
        %Pick out left and right profiles
        Ne_1d_left_temp = fliplr(den_int(1:centroid_ind_tmp));
        Ne_1d_right_temp = den_int(centroid_ind_tmp+1:end);
        [Ne_min_l,ind_Ne_min_l] = min(Ne_1d_left_temp);
        [Ne_min_r,ind_Ne_min_r] = min(Ne_1d_right_temp);
        
        %Truncate profiles if the minimum density occurs before the edge of
        %the hologram
        
        if ind_Ne_min_l ~= length(Ne_1d_left_temp)
            truncate_l(jcol,jcen) = length(Ne_1d_left_temp) - ind_Ne_min_l;
            Ne_1d_left{jcol,jcen} = Ne_1d_left_temp(1:ind_Ne_min_l);
        else
            truncate_l(jcol,jcen) = 0;
            Ne_1d_left{jcol,jcen} = Ne_1d_left_temp;
        end
        
        if ind_Ne_min_r ~= length(Ne_1d_right_temp)
            truncate_r(jcol,jcen) = length(Ne_1d_right_temp) - ind_Ne_min_r;
            Ne_1d_right{jcol,jcen} = Ne_1d_right_temp(1:ind_Ne_min_r);
        else
            truncate_r(jcol,jcen) = 0;
            Ne_1d_right{jcol,jcen} = Ne_1d_right_temp;
        end
        
        %abel invert
        ne_1d_left{jcol,jcen} = abel_invert(y_twin_red_adj,[Ne_1d_left{jcol,jcen} -...
            min(Ne_1d_left{jcol,jcen})]');
        ne_1d_right{jcol,jcen} = abel_invert(y_twin_red_adj,[Ne_1d_right{jcol,jcen}-...
            min(Ne_1d_right{jcol,jcen})]');
        
        centroid_ind_track(jcol,jcen) = centroid_ind_tmp;
        centroid_ind_tmp = centroid_ind_tmp + 1;
        
        %find the length of each density vector
        len_l = length(ne_1d_left{jcol,jcen});
        len_r=length(ne_1d_right{jcol,jcen});
        
        %Truncate the longer side to compare the L2 norm
        if len_l > len_r
            ne_1d_left_compare{jcol,jcen} = ne_1d_left{jcol,jcen}(1:len_r);
            ne_1d_right_compare{jcol,jcen} = ne_1d_right{jcol,jcen};
        elseif len_r > len_l
            ne_1d_right_compare{jcol,jcen} = ne_1d_right{jcol,jcen}(1:len_l);
            ne_1d_left_compare{jcol,jcen} = ne_1d_left{jcol,jcen};
        else
            ne_1d_right_compare{jcol,jcen} = ne_1d_right{jcol,jcen};
            ne_1d_left_compare{jcol,jcen} = ne_1d_left{jcol,jcen};
        end
        
        n_val = 2E23;   %value to add to the shorter profile
        ind_edge = -n_tries/2; %index to multipy n_val
        
        %add the specified density value to the shorter profile
        for jd = 1:n_tries
            if len_l>len_r
                ne_1d_left_add{jcol,jcen,jd} = ne_1d_left_compare{jcol,jcen};
                ne_1d_right_add{jcol,jcen,jd} = ne_1d_right_compare{jcol,jcen}(:) + n_val*ind_edge;
            else
                ne_1d_right_add{jcol,jcen,jd} = ne_1d_right_compare{jcol,jcen};
                ne_1d_left_add{jcol,jcen,jd} = ne_1d_left_compare{jcol,jcen}(:) + n_val*ind_edge;
            end
            
            ind_edge = ind_edge + 1; %increment density multiplier
            
            n_ind_track(jd) = n_val*ind_edge;
            
            c_dens_diff{jcol,jcen,jd} = abs(ne_1d_left_add{jcol,jcen,jd}(1:5) - ne_1d_right_add{jcol,jcen,jd}(1:5));
            c_norm_dens(jcol,jcen,jd) = norm(c_dens_diff{jcol,jcen,jd},2)/5;
        end
        
        %Find the minimum of the L2 norm to only take one density profile
        %per centroid. This tells you how much density to add/subtract from
        %each profile.
        [cmin(jcol,jcen),cind(jcol,jcen)] = min(c_norm_dens(jcol,jcen,:));
        
    end
    
    %Find the minimum of the L2 norm for each profile to select the optimal
    %centroid
    [f_min(jcol),f_ind(jcol)] = min(cmin(jcol,:));
    
end