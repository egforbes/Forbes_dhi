%This code is designed to be used as a diagnostic
indf = [];
for jlook = 200;
   den_int = den_int_full(jlook,:); 
   [Ne_max(jlook),ind_Ne_max(jlook)] = max(den_int);
   
    centroid_ind_tmp = ind_Ne_max(k);
    
    Ne_1d_left{k} = fliplr(den_int(1:centroid_ind_tmp));
    Ne_1d_right{k} = den_int(centroid_ind_tmp+1:end);
end


