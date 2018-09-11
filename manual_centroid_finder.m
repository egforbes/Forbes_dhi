clear all; close all; clc;
% load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\',num2str(171005),'\Ne_',num2str(171005023),'_long.mat'))
% load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\',num2str(171005),'\Ne_',num2str(171005023),'.mat'))
% load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\backup_good','\Ne_',num2str(171005023),'_backup.mat'))
load(strcat('W:\Users\Eleanor_Forbes\Reconstructions\',num2str(171005),'\Ne_',num2str(171005023),'sm.mat'))
% den_int = den_int_full(45,:); %If using the full reconstruction, need den_int_full(:,x)
den_int_full = den_int_full(2:end,:);
den_int = den_int_full(:,5).';
[Ne_max, ind_Ne_max] = max(den_int);  %Find the location of the maximum
num_of_centroids = 6; %Specify number of centroid guesses

n_tries = 16; %this is how many times you add density to the short side of the profile. Must be EVEN or it will barf

%%
centroid_ind_tmp = ind_Ne_max-num_of_centroids; %Starting location of the centroid
for j = 1:2*num_of_centroids
    
    Ne_1d_left_temp = fliplr(den_int(1:centroid_ind_tmp));
% Ne_1d_left_temp = flipud(den_int(1:centroid_ind_tmp));
    Ne_1d_right_temp = den_int(centroid_ind_tmp+1:end);
    [Ne_min_l,ind_Ne_min_l] = min(Ne_1d_left_temp);
    [Ne_min_r,ind_Ne_min_r] = min(Ne_1d_right_temp);
    
    
    
    if ind_Ne_min_l ~= length(Ne_1d_left_temp)
        truncate_l(j) = length(Ne_1d_left_temp) - ind_Ne_min_l;
        Ne_1d_left{j} = Ne_1d_left_temp(1:ind_Ne_min_l);
    else
        truncate_l(j) = 0;
        Ne_1d_left{j} = Ne_1d_left_temp;
    end
    
    if ind_Ne_min_r ~= length(Ne_1d_right_temp)
        truncate_r(j) = length(Ne_1d_right_temp) - ind_Ne_min_r;
        Ne_1d_right{j} = Ne_1d_right_temp(1:ind_Ne_min_r);
    else
        truncate_r(j) = 0;
        Ne_1d_right{j} = Ne_1d_right_temp;
    end
    

    ne_1d_left{j} = abel_invert(y_twin_red_adj,[Ne_1d_left{j}-...
        min(Ne_1d_left{j})]');
    ne_1d_right{j} = abel_invert(y_twin_red_adj,[Ne_1d_right{j}-...
        min(Ne_1d_right{j})]');
    
    centroid_ind_track(j) = centroid_ind_tmp;
    centroid_ind_tmp = centroid_ind_tmp + 1;
    
    len_l= length(ne_1d_left{j});
    len_r=length(ne_1d_right{j});
    
        figure(1)
    subplot(4,5,j)
    plot([1:centroid_ind_track(j)],fliplr(Ne_1d_left{j}),'ko','Linewidth',2)
    hold on
    plot([(centroid_ind_track(j)+1): len_l+len_r],Ne_1d_right{j},'bs','Linewidth',2)
    hold off
    
    figure(2)
    subplot(4,5,j)
    plot([1:centroid_ind_track(j)],flipud(ne_1d_left{j}),'Linewidth',2)
    hold on
    plot([(centroid_ind_track(j)+1): len_l+len_r], ne_1d_right{j},'Linewidth',2)
    hold off
    
    
    %Truncate the longer side to compare the L2 norm
    
    
    if len_l > len_r
        ne_1d_left_compare{j} = ne_1d_left{j}(1:len_r);
        ne_1d_right_compare{j} = ne_1d_right{j};
    elseif len_r > len_l
        ne_1d_right_compare{j} = ne_1d_right{j}(1:len_l);
        ne_1d_left_compare{j} = ne_1d_left{j};
    else
        ne_1d_right_compare{j} = ne_1d_right{j};
        ne_1d_left_compare{j} = ne_1d_left{j};
    end
    
    n_val = 2E23;

ind_edge = -n_tries/2;
    fignum = j + 2;
    
    for jdiff = 1:n_tries
        if len_l>len_r
            ne_1d_left_add{j,jdiff} = ne_1d_left_compare{j};
            ne_1d_right_add{j,jdiff} = ne_1d_right_compare{j}(:) + n_val*ind_edge;
        else
            ne_1d_right_add{j,jdiff} = ne_1d_right_compare{j};
            ne_1d_left_add{j,jdiff} = ne_1d_left_compare{j}(:) + n_val*ind_edge;
        end
        ind_edge = ind_edge + 1;
%         figure(fignum)
%         title(num2str(j))
%         subplot(2,5,jdiff)
%         plot(ne_1d_left_add{j,jdiff},'Linewidth',2)
%         hold on
%         plot(ne_1d_right_add{j,jdiff},'Linewidth',2)
%         hold off
        n_ind_track(jdiff) = n_val*ind_edge;
        
        mini_dens_diff{j,jdiff} = abs(ne_1d_left_add{j,jdiff}-ne_1d_right_add{j,jdiff});
        mini_norm_dens(j,jdiff) = norm(mini_dens_diff{j,jdiff},2)/length(mini_dens_diff{j,jdiff});
        
        fmin_dens_diff{j,jdiff} = abs(ne_1d_left_add{j,jdiff}(1:7)-ne_1d_right_add{j,jdiff}(1:7));
        fmin_norm_dens(j,jdiff) = norm(fmin_dens_diff{j,jdiff},2)/length(mini_dens_diff{j,jdiff}(1:5));
        
        mmin_dens_diff(j,jdiff) = abs(ne_1d_left_add{j,jdiff}(1)-ne_1d_right_add{j,jdiff}(1));
%         mmin_norm_dens(j,jdiff) = norm(mmin_dens_diff{j,jdiff},2);
        
    end
    
    
    
    dens_diff{j} = abs(ne_1d_left_compare{j} - ne_1d_right_compare{j});
    norm_dens(j) = norm(dens_diff{j},2)/(length(dens_diff{j}));
%     
%     
%     %now compare the L2 norm for only the central part (first 15 points)
    dens_diff_t{j} = abs(ne_1d_left_compare{j}(1:15) - ne_1d_right_compare{j}(1:15));
    norm_dens_t(j) = norm(dens_diff_t{j},2)/(length(dens_diff{j}));
%     
%     
%     
%     
%     %First five points
    dens_diff_s{j} = abs(ne_1d_left_compare{j}(1:7) - ne_1d_right_compare{j}(1:7));
    norm_dens_s(j) = norm(dens_diff_s{j},2)/(length(dens_diff{j}));
    
    dens_diff_sm{j} = abs(ne_1d_left_compare{j}(1) - ne_1d_right_compare{j}(1));
    norm_dens_sm(j) = norm(dens_diff_t{j},2)/(length(dens_diff{j}));
    
end
% 
% figure(34)
% plot(n_ind_track,mini_norm_dens(6,:),'ko')
% hold on
% plot(n_ind_track,fmin_norm_dens(8,:),'bs')
% hold off
%%
for jl = 1:size(mini_norm_dens,1)
   [mini(jl),indmini(jl)] = min(mini_norm_dens(jl,:)); 
   [fmin(jl),findmin(jl)] = min(fmin_norm_dens(jl,:));
   [emin(jl),eind(jl)] = min(mmin_dens_diff(jl,:));
end



[exmin, indexmin] = min(mini);
[fex, fexmin] = min(fmin);
[bex,bexmin] = min(emin);
figure(666)
plot(centroid_ind_track,norm_dens,'ro','Linewidth',2)
hold on
plot(centroid_ind_track,norm_dens_t,'bs','Linewidth',2)
plot(centroid_ind_track,norm_dens_s,'k.','Linewidth',2)
% plot(centroid_ind_track,norm_dens_sm,'m*','Linewidth',2)
xlabel('Centroid Choice (Pixel)')
ylabel('Normalized Density Difference')
title('Manual Centroid Finder Results')
set(gca,'Fontsize',14)
legend('All Points','15 Points','5 points')

figure()
subplot(1,3,1)
plot(ne_1d_left_add{fexmin,findmin(fexmin)},'bs','linewidth',2)
hold on
plot(ne_1d_right_add{fexmin,findmin(fexmin)},'ko','linewidth',2)
hold off
legend('left','right')
set(gca,'Fontsize',14)
title('Center Matching')
subplot(1,3,2)
plot(ne_1d_left_add{indexmin,indmini(indexmin)},'bs','linewidth',2)
hold on
plot(ne_1d_right_add{indexmin,indmini(indexmin)},'ko','linewidth',2)
hold off
legend('left','right')
set(gca,'Fontsize',14)
title('Full Profile Matching')
subplot(1,3,3)
plot(ne_1d_left_add{bexmin,eind(bexmin)},'bs','linewidth',2)
hold on
plot(ne_1d_right_add{bexmin,eind(bexmin)},'ko','linewidth',2)
hold off
legend('left','right')
set(gca,'Fontsize',14)
title('Center Point Matching')


%%
[mind, indmin] = min(norm_dens);
compvec = [flipud(ne_1d_left{indmin}); ne_1d_right{indmin}];

[smin,sindmin] = min(norm_dens_sm);
ccvec = [flipud(ne_1d_left{sindmin}); ne_1d_right{sindmin}];

figure(56)
subplot(1,2,1)
plot([1:centroid_ind_track(indmin)],flipud(ne_1d_left{indmin}),'bs','Linewidth',2)
hold on
plot([centroid_ind_track(indmin) + 1: length(compvec)],ne_1d_right{indmin},'ro','Linewidth',2)
hold off
legend('Left Profile','Right Profile')
title('All Points "Best" Profile')

subplot(1,2,2)
plot([1:centroid_ind_track(sindmin)],flipud(ne_1d_left{sindmin}),'bs','Linewidth',2)
hold on
plot([centroid_ind_track(sindmin) + 1: length(compvec)],ne_1d_right{sindmin},'ro','Linewidth',2)
hold off
legend('Left Profile','Right Profile')
title('Center Matching "Best" Profile')

%%
%Now we need to put it all together, add the density value from the center
%matching