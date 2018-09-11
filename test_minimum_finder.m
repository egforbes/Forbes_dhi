load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\',num2str(171005),'\Ne_',num2str(171005023),'.mat'))
den_int_full = den_int_full(1:75,:);
num_of_centroids = 2;
num_of_n_edge = 1;

 for k = 1:size(den_int_full,1)
% for k = 65

    den_int = den_int_full(k,:);
    [Ne_max(k),ind_Ne_max(k)] = max(den_int);
    
    centroid_ind_tmp = ind_Ne_max(k)-num_of_centroids;
    if centroid_ind_tmp <=0
        den_num_l{k} = NaN(size(den_int_full,1),1);%[];
        den_num_r{k} = NaN(size(den_int_full,1),1);
        centroid_abs(k) = 1;
        centroid(k) = 0;
        % for centroids at the upper boundary:
    elseif centroid_ind_tmp >= size(den_int_full,1)
        den_num_l{k} = NaN(size(den_int_full,1),1);
        den_num_r{k} = NaN(size(den_int_full,1),1);%[];clc
        centroid_abs(k) = size(den_int_full,1);
        centroid(k) = 0;
    else
        
        
        % for centroids near (within 2*num_of_centroids of) the upper
        % boundary:
        if centroid_ind_tmp >= size(den_int_full,1)-2*num_of_centroids
            centroids_iterations = 2*num_of_centroids-...
                (centroid_ind_tmp-(size(den_int_full,1)-2*...
                num_of_centroids));
            % for centroids near the lower boundary:
        else
            centroids_iterations = 2*num_of_centroids;
        end
    end
   
    
    for j = 1:centroids_iterations
        clear Ne_1d_left_temp Ne_1d_right_temp
        Ne_1d_left_temp = fliplr(den_int(1:centroid_ind_tmp));
        Ne_1d_right_temp = den_int(centroid_ind_tmp+1:end);
        [Ne_min_l,ind_Ne_min_l] = min(Ne_1d_left_temp);
        [Ne_min_r,ind_Ne_min_r] = min(Ne_1d_right_temp);
        
        
%                     figure(1)
%                     plot(Ne_1d_right_temp)
%                     hold on
%                     plot(ind_Ne_min_r,Ne_min_r,'ko')
%                     hold off
%                     pause(0.1)
        
        if ind_Ne_min_l ~= length(Ne_1d_left_temp)
            truncate_l(k,j) = length(Ne_1d_left_temp) - ind_Ne_min_l;
            Ne_1d_left{k,j} = Ne_1d_left_temp(1:ind_Ne_min_l);
        else
            truncate_l(k,j) = 0;
            Ne_1d_left{k,j} = Ne_1d_left_temp;
        end
        
        if ind_Ne_min_r ~= length(Ne_1d_right_temp)
            truncate_r(k,j) = length(Ne_1d_right_temp) - ind_Ne_min_r;
            Ne_1d_right{k,j} = Ne_1d_right_temp(1:ind_Ne_min_r);
        else
            truncate_r(k,j) = 0;
            Ne_1d_right{k,j} = Ne_1d_right_temp;
        end
        ne_1d_left{k,j} = abel_invert(y_twin_red_adj,[Ne_1d_left{k,j}-...
        min(Ne_1d_left{k,j})]');
        ne_1d_right{k,j} = abel_invert(y_twin_red_adj,[Ne_1d_right{k,j}-...
        min(Ne_1d_right{k,j})]');
%     figure(1)
%     plot(Ne_1d_right_temp)
%     hold on
%     plot(ind_Ne_min_r,Ne_min_r,'ko')
%     plot(Ne_1d_right{k,j},'.-')
%     hold off
%     pause(0.1)
%     

% Iterate the index to prepare the code to move to the next
                % centroid value on the next loop iteration:
                centroid_ind_track(k,j) = centroid_ind_tmp;
                centroid_ind_tmp = centroid_ind_tmp + 1;


    
                    % When the centroid does not divide the Ne profile in half,
            % the longer of the two sides must be truncated to match
            % the length of the shorter side.  The following case
            % statements conduct this truncation:
            len_left(k,j) = length(Ne_1d_left{k,j});
            len_right(k,j) = length(Ne_1d_right{k,j});
            
            if len_left(k,j)>len_right(k,j)
                ne_1d_left_trunc{k,j} = ne_1d_left{k,j}(1:...
                    len_right(k,j));
                ne_1d_right_trunc{k,j} = ne_1d_right{k,j};
            elseif len_left(k,j)<len_right(k,j)
                ne_1d_left_trunc{k,j} = ne_1d_left{k,j};
                ne_1d_right_trunc{k,j} = ne_1d_right{k,j}(1:...
                    len_left(k,j));
            else
                ne_1d_left_trunc{k,j} = ne_1d_left{k,j};
                ne_1d_right_trunc{k,j} = ne_1d_right{k,j};
            end
            
            % Identify the nominal n_edge_0 value as the density in the
            % longer profile at the radius of the edge of the short
            % profile:
            if len_left(k,j)>len_right(k,j) && len_right(k,j)~=0
                n_edge_0(k,j) = ne_1d_left_trunc{k,j}(len_right(k,j));
            elseif len_left(k,j)<len_right(k,j) && len_left(k,j)~=0
                n_edge_0(k,j) = ne_1d_right_trunc{k,j}(len_left(k,j));
            else
                n_edge_0(k,j) = 0;
            end
            
            % Initialize n_edge_tmp to enable determination of
            % centroid_abs (the centroid locations in absolute
            % pixel-space):
            n_edge_ind_tmp = -num_of_n_edge;
            num_of_n_edge = 2;
            n_edge_delta = 5e21;
              % Iterate through n_edge values:
            for jj = 1:2*num_of_n_edge
                % populating list of all n_edge attempted:
                n_edge_val(jj) = n_edge_ind_tmp*n_edge_delta;
                
                % Iterating through n_edge values for the shorter
                % profile:
                if len_left(k,j)>len_right(k,j)
                    ne_1d_left_compare{k,j,jj} = ne_1d_left_trunc{k,j};
                    ne_1d_right_compare{k,j,jj} = ...
                        ne_1d_right_trunc{k,j}+...
                        n_edge_delta*n_edge_ind_tmp;
                elseif len_left(k,j)<len_right(k,j)
                    ne_1d_left_compare{k,j,jj} = ...
                        ne_1d_left_trunc{k,j}+...
                        n_edge_delta*n_edge_ind_tmp;
                    ne_1d_right_compare{k,j,jj} = ...
                        ne_1d_right_trunc{k,j};
                else
                    if ne_1d_left_trunc{k,j}(end) > ...
                            ne_1d_right_trunc{k,j}(end)
                        ne_1d_left_compare{k,j,jj} = ...
                            ne_1d_left_trunc{k,j}+...
                            n_edge_delta*n_edge_ind_tmp;
                        ne_1d_right_compare{k,j,jj} = ...
                            ne_1d_right_trunc{k,j};
                    else
                        ne_1d_left_compare{k,j,jj} = ...
                            ne_1d_left_trunc{k,j};
                        ne_1d_right_compare{k,j,jj} = ...
                            ne_1d_right_trunc{k,j}+...
                            n_edge_delta*n_edge_ind_tmp;
                    end
                end
                
                % Iterate the index to prepare the code to move to the
                % next n_edge value on the next loop iteration:
                n_edge_ind_track(k,j,jj) = n_edge_ind_tmp;
                n_edge_ind_tmp = n_edge_ind_tmp + 1;
                
                % Computing the L2 norm of the difference in the sides
                % of the ne profile:
                dens_diff{k,j,jj} = abs(ne_1d_left_compare{k,j,jj}-...
                    ne_1d_right_compare{k,j,jj});
                norm_dens(k,j,jj) = norm(dens_diff{k,j,jj},2)/...
                    (length(dens_diff{k,j,jj}));
            end
            
            % Changing all zero-valued elements in norm_dens to NaNs
            % (allows the code to ignore them when searching for the
            % minimum)
            ind_zeros = find(norm_dens(k,j,:)==0);
            norm_dens(k,j,ind_zeros) = NaN;
            
            % Find the minimum L2 norm value for each combination of
            % cross-section and centroid.  (n_edge_intermediate records
            % the index location of n_edge yielding minimum L2 norm)
            [norm_intermediate(k,j),n_edge_intermediate(k,j)] = ...
                min(norm_dens(k,j,:),[],3);
         end
        
        % Find the minimum L2 norm value for each cross-section.
        % (centroid records the relative index location of the centroid
        % location yielding minimum L2 norm)
        [norm_min_val(k),centroid(k)] = min(norm_intermediate(k,:)...
            ,[],2);
        n_edge(k) = n_edge_intermediate(k,centroid(k));
        
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
            if ne_1d_left_trunc{k,j}(end) > ne_1d_right_trunc{k,j}(end)
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
        end
        
        
        % Computing the centroid indices in absolute pixel coordinates:
        if centroid_ind_track(k,centroid(k)) == 0
            centroid_abs(k) = size(den_int_full,1);
        else
            centroid_abs(k) = centroid_ind_track(k,centroid(k));
        end
    if truncate_l(k,centroid(k)) ~= 0
        den_num_l_fin{k} = [den_num_l{k}; NaN(1,truncate_l(k,centroid(k)))];
    end
    if truncate_r(k,centroid(k)) ~= 0;
        den_num_r_fin{k} = [den_num_r{k}; NaN(truncate_r(k,centroid(k)),1)];
    end
   
    
 end   

 %%
% Compiling the inverted number density data to form the
% den_num_full matrix:
    den_num_full = zeros(size(den_int_full));
    num_of_cross_sect = size(den_int_full,1);

    for k = 1:num_of_cross_sect
        if centroid_abs(k) == 1 || centroid_abs(k) == size(den_int_full,1)
            den_num_full(k,:) = zeros(size(den_int_full,1),1);
            
        else
            den_num_full(k,(1+centroid_abs(k)-length(den_num_l{k})):...
                centroid_abs(k)) = flipud(den_num_l{k});
            den_num_full(k,(1+centroid_abs(k)):(centroid_abs(k)+...
                length(den_num_r{k}))) = den_num_r{k};
        end
    end
    
  pcolor(x_twin_red_adj(1:75),y_twin_red_adj, den_num_full.')
  shading interp
  colorbar
  
   
