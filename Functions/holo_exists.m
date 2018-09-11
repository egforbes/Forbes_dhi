% holo_exists:  checks to see if the holograms are not blanks.

function [skip] = holo_exists(shotnum_base,shotnum_def)

%%% Starting with ZaP's shot number convention:
SN_def = num2str(shotnum_def);
SN_base = num2str(shotnum_base);

%%% Converting to Canon Camera shot number convention:
year_def = num2str(2000+str2num(SN_def(1:2)));
month_def = num2str(SN_def(3:4));
day_def = num2str(SN_def(5:6));

year_base = num2str(2000+str2num(SN_base(1:2)));
month_base = num2str(SN_base(3:4));
day_base = num2str(SN_base(5:6));

date_def = [year_def,'_', month_def,'_', day_def];
date_base = [year_base,'_',month_base,'_',day_base];
date_plot = [year_def, month_def, day_def];

date_plot = 999;
shot_def = 999;
shot_base = 999;

shot_def = SN_def(7:end);
shot_base = SN_base(7:end);

hol_def_rgb = imread(['H:\Raw Holograms\',date_def,'\',...
    date_def, '_', shot_def,'.jpg']);
hol_base_rgb = imread(['H:\Raw Holograms\',date_base,'\',...
    date_base, '_', shot_base,'.jpg']);

%%% Converting from RGB to grayscale (rgb2gray is not a supported command
%%% without the image processing toolbox):
rgb_red = 0.3;
rgb_green = 0.59;
rgb_blue = 0.11;

hol_def = double(hol_def_rgb(:,:,1).*rgb_red + hol_def_rgb(:,:,2).*rgb_green + hol_def_rgb(:,:,3).*rgb_blue)';
hol_base = double(hol_base_rgb(:,:,1).*rgb_red + hol_base_rgb(:,:,2).*rgb_green + hol_base_rgb(:,:,3).*rgb_blue)';

hol_def_int = mean(mean(hol_def));
hol_base_int = mean(mean(hol_base));

% check to make sure each hologram actually exists (is not a blank):

if hol_def_int<20|hol_base_int<20
    skip = 1;
else
    skip = 0;
end