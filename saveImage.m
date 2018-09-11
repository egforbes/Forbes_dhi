fileName = '../../data/Plasma_Shot_28.JPG';
%fileName = '../../data/Base_Shot_26.JPG';
hol_base_rgb = imread(fileName);

rgb_red = 0.3;
rgb_green = 0.59;
rgb_blue = 0.11;
    
hol_base = double(hol_base_rgb(:,:,1).*rgb_red + hol_base_rgb(:,:,2).*rgb_green + hol_base_rgb(:,:,3).*rgb_blue)';

newFileName = strrep(fileName,'JPG','dat');

fidR = fopen(newFileName,'w');
fwrite(fidR,hol_base,'double');
fclose(fidR);

image(hol_base)
