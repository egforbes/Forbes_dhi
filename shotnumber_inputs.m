% shotnumber_inputs.m
% note that the final shot in each group needs to be a deformed shot
% for the deformed-baseline pairing of shot_base_pair.m to work.

function [shots,d,xmin,xmax,ymin,ymax,sign_twin] = shotnumber_inputs(import,varargin)

% using varargin to accept input identifying if the user desires
% deformed-baseline or baseline-baseline pairs
switch(nargin),
    case 0,
        whichPair = 1; % if second input is not provided, the program seeks
        %deformed-baseline pairs by default
    case 1,
        whichPair = 1; % if second input is not provided, the program seeks
        %deformed-baseline pairs by default
    case 2,
        whichPair = varargin{1}; % a second input value of 1 tells the
        %program to seek for deformed-baseline pairs, a value of 0 tells
        %the program to seek for baseline-baseline pairs
end

switch(import) % this switch statement gives the user some flexibility to
    % toggle between different groups of shots for analysis.  Usually, case
    % 2 is used.
    case 1
        
        shots_tmp = shot_base_pair(shotnum,whichPair);
        d(1:size(shots_tmp,1)) = 0.5;
        xmin(1:size(shots_tmp,1)) = 0.0446;
        xmax(1:size(shots_tmp,1)) = 0.0555;
        ymin(1:size(shots_tmp,1)) = 0.0237;
        ymax(1:size(shots_tmp,1)) = 0.0429;
        sign_twin(1:size(shots_tmp,1)) = -1;
        shots = shots_tmp;
        
    case 2
        
        shots = 171005000 + [21 23];
% %         shots = 170927000 + [26 32];
% %         shots = 170927000 + [26 30];
        shots_tmp=shot_base_pair(shots,whichPair);
        d=0.45;
        xmin=0.04259;
        xmax=0.05255;%xmax=0.05355;
        ymin=0.02632;
        ymax=0.04196;
        sign_twin=-1;
        shots=shots_tmp;
%         
% shots = 180215000 + [12 13];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.55;
%         xmin=0.04864;  
%         xmax=0.06363;   
%         ymin=0.03129; 
%         ymax=0.05413; 
%         sign_twin=-1;
%         shots=shots_tmp;

        
%         shots = 171005000+[7:23];
% shots = 171005000 + [21 23];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.43;
%         xmin=0.03862;  %0.03772;
%         xmax=0.05061;   %0.05058;
%         ymin=0.02551; %0.02528;
%         ymax=0.03824; %0.03804;
%         sign_twin=-1;
%         shots=shots_tmp;

%          shots = 171005000+[7:24];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.55;
%         xmin=0.04923;
%         xmax=0.06407;
%         ymin=0.03042;
%         ymax=0.04879;
%         sign_twin=-1;
%         shots=shots_tmp;

%         shots = 170927000+[37 39];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.6;
%         xmin=0.05476;
%         xmax=0.06942;
%         ymin=0.03094;
%         ymax=0.05118;
%         sign_twin=-1;
%         shots=shots_tmp;
%         
%         shots = 170901000+[7:45];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.55;
%         xmin=0.05016;
%         xmax=0.06496;
%         ymin=0.0282;
%         ymax=0.05084;
%         sign_twin=-1;
%         shots=shots_tmp;
        
%         shots = 170901000+[7:45];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.5;
%         xmin=0.04639;
%         xmax=0.06159;
%         ymin=0.02436;
%         ymax=0.04806;
%         sign_twin=-1;
%         shots=shots_tmp;
        
%           shots = 170824000+[7:53];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.5;
%         xmin=0.04615;
%         xmax=0.0614;
%         ymin=0.02478;
%         ymax=0.04792;
%         sign_twin=-1;
%         shots=shots_tmp;
%         shots = 170504000+[42:63]; 
%         shots = 170510000+[1:9];


%         shots = 170622000+[8 27];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.65;
%       xmin=0.0571;
%         xmax=0.07264;
%         ymin=0.03483;
%         ymax=0.05797;
%         sign_twin=-1;
%         shots=shots_tmp;


%         shots = 170608000+[11:24];
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.65;
%         xmin=0.0571;
%         xmax=0.07264;
%         ymin=0.03483;
%         ymax=0.05797;
%         sign_twin=-1;
%         shots=shots_tmp;
        
%         shots = 170420000+[26:28]; 
%         shots_tmp=shot_base_pair(shots,whichPair);
%         d=0.6;
%         xmin=0.05223;
%         xmax=0.06757;
%         ymin=0.03155;
%         ymax=0.05258;
%  d = 0.5;
%                 xmin = 0.0434;
%                 xmax = 0.0583;
%                 ymin = 0.0245;
%                 ymax = 0.0440;
%         sign_twin=-1;
%         shots=shots_tmp;
        
        %         shots = 160823000+[15,16]; %shots = [160803000+[5:28]];
        %         shots_tmp = shot_base_pair(shots,whichPair);
        %         d = 0.5;
        %         xmin = 0.0438;
        %         xmax = 0.0565;
        %         ymin = 0.0245;
        %         ymax = 0.0450;
        % %         xmin = 0.044;
        % %         xmax = 0.0565;
        % %         ymin = 0.0233;
        % %         ymax = 0.0445;
        %         sign_twin = -1;
        %         shots = shots_tmp;
        
        %     shotnum = 160323000+[7:50];
        %     shots = shot_base_pair(shotnum,whichPair);
        %     d = 1.2;
        %     xmin = 0.092;
        %     xmax = 0.105;
        %     ymin = 0.0737;
        %     ymax = 0.0943;
        %     sign_twin = -1;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         shots = [160127047 160127045];%+[8:70];
        %         d = 0.88;
        %         xmin = 0.0707;
        %         xmax = 0.0827;
        %         ymin = 0.0570;
        %         ymax = 0.0761;
        %         sign_twin = -1;
        
        %         shots = [150915008 150915007];
        %         d =.40; %      % Reconstruction distance [m]
        %         xmin = 0.0058;  %0.006;
        %         xmax = 0.0188;  %0.019;
        %         ymin = 0.01685;  %0.015;
        %         ymax = 0.032;  %0.035;
        %         sign_twin = -1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         shots = [160622000+[6:26]];
        %         shots = [160524000+[20,21]];
        %         shots_tmp = shot_base_pair(shots,whichPair);
        %         d = 0.5;
        % %         xmin = 0.0434;
        % %         xmax = 0.0583;
        % %         ymin = 0.0245;
        % %         ymax = 0.0440;
        %         xmin = 0.044;
        %         xmax = 0.0565;
        %         ymin = 0.0233;
        %         ymax = 0.0445;
        %         sign_twin = -1;
        %         shots = shots_tmp;
        
        
        % RSI 2016 Paper shots for analysis:
%                 shots = [160524000+[20,21]];%[21,20]];
        %Forbes analysis for 170420
% %                 shots = [170420000+[8:34]];
% %                 shots_tmp = shot_base_pair(shots,whichPair);
% %                 d = 0.5;
% %                 xmin = 0.0434;
% %                 xmax = 0.0583;
% %                 ymin = 0.0245;
% %                 ymax = 0.0440;


%                 xmin = 0.0446;
%                 xmax = 0.0555;
%                 ymin = 0.0237;
%                 ymax = 0.0429;


%                 sign_twin = -1;
%                 shots = shots_tmp;
%         
        %         shots = [160524018];%[21,20]];
        % %         shots_tmp = shot_base_pair(shots,whichPair);
        %         d = 0.5;
        % %         xmin = 0.0434;
        % %         xmax = 0.0583;
        % %         ymin = 0.0245;
        % %         ymax = 0.0440;
        %         xmin = 0.0446;
        %         xmax = 0.0555;
        %         ymin = 0.0237;
        %         ymax = 0.0429;
        %         sign_twin = -1;
        % %         shots = shots_tmp;
        
        %         shots = 160823016;%+[5:45];
        % %         shots = shot_base_pair(shotnum,whichPair);
        %         d = 1.2;
        %         xmin = 0.092;
        %         xmax = 0.105;
        %         ymin = 0.0737;
        %         ymax = 0.0943;
        %         sign_twin = -1;
        
        
%         shots = 160310013;%+[5:45];
%         %         shots = shot_base_pair(shotnum,whichPair);
%         d = 1.2;
%         xmin = 0.092;
%         xmax = 0.105;
%         ymin = 0.0737;
%         ymax = 0.0943;
%         sign_twin = -1;
        
        
end




% if import==1

%     shotnum = 160323000+[7:50];
%     shots = shot_base_pair(shotnum,whichPair);
%     d = 1.2;
%     xmin = 0.092;
%     xmax = 0.105;
%     ymin = 0.0737;
%     ymax = 0.0943;
%     sign_twin = -1;

%     shotnum = 160316000+[12:45];
%     shots = shot_base_pair(shotnum);
%     d = 1.2;
%     xmin = 0.092;
%     xmax = 0.105;
%     ymin = 0.0737;
%     ymax = 0.0943;
%     sign_twin = -1;

%     shotnum = 160310000+[24:25];
%     shots = shot_base_pair(shotnum);
%     d = 0.85;
%     xmin = 0.0671;
%     xmax = 0.0799;
%     ymin = 0.0538;
%     ymax = 0.0740;
%     sign_twin = -1;


%     shotnum = 160224000+[7:55];
%     shots = shot_base_pair(shotnum);
%     d = 0.85;
%     xmin = 0.0671;
%     xmax = 0.0799;
%     ymin = 0.0538;
%     ymax = 0.0740;
%     sign_twin = -1;

% %     shotnum = 160203000+[3:41];
%     shots = [160203037 160203036];
% %     shots = shot_base_pair(shotnum);
%     d = 0.8128;
%     xmin = 0.0639;
%     xmax = 0.0769;
%     ymin = 0.0508;
%     ymax = 0.0711;
%     sign_twin = -1;

% %     shotnum = 160202000+[5:30];
%     shots = [160202007 160202005];
% %     shots = shot_base_pair(shotnum);
%     d = 0.8128;
%     xmin = 0.0639;
%     xmax = 0.0769;
%     ymin = 0.0508;
%     ymax = 0.0711;
%     sign_twin = -1;

%     shotnum = 151028000+[9:77];
%
%     shots = shot_base_pair(shotnum);
%     d = .38;
%     xmin = 0.0345;
%     xmax = 0.0458;
%     ymin = 0.0209;
%     ymax = 0.0403;
%     sign_twin = -1;



%     shots = [
%     151029008 151029007;
%     151029011 151029010;
%     151029012 151029010;
%     151029013 151029010;
%     151029015 151029014;
%     151029017 151029016;
%     151029018 151029016;
%     151029020 151029019;
%     151029022 151029021;
%     151029023 151029021;
%     151029025 151029024;
%     151029027 151029026;
%     151029028 151029026;
%     151029030 151029029;
%     151029032 151029031;
%     151029033 151029031;
%     151029035 151029034;
%     151029037 151029036;
% ];
%     d = .38;
%     xmin = 0.0340;
%     xmax = 0.0458;
%     ymin = 0.0209;
%     ymax = 0.0403;
%     sign_twin = -1;


%     shots = [
%     151030007 151030006;
%     151030008 151030006;
%     151030009 151030006;
%     151030011 151030010;
%     151030012 151030010;
%     151030013 151030010;
%     151030014 151030010;
%     151030015 151030010;
%     151030017 151030016;
%     151030018 151030016;
%     151030020 151030019;
%     151030021 151030019;
%     151030023 151030022;
%     151030024 151030022;
%     151030025 151030022;
%     151030027 151030026;
%     151030028 151030026;
%     151030029 151030026;
%     151030031 151030030;
%     151030032 151030030;
%     151030033 151030030;
% ];
%     d = .38;
%     xmin = 0.0340;
%     xmax = 0.0458;
%     ymin = 0.0209;
%     ymax = 0.0403;
%     sign_twin = -1;




% shots = [
%     160120011 160120010;
%     160120015 160120014;
%     160120017 160120016;
%     160120019 160120018;
%     160120021 160120020;
%     160120023 160120022;
%     160120025 160120024;
%     160120027 160120026;
%     160120029 160120028;
%     160120031 160120030;
%     160120033 160120032;
%     160120035 160120034;
%     160120037 160120036;
% ];
%     d = .88;
%     xmin = 0.0707;
%     xmax = 0.0838;
%     ymin = 0.0570;
%     ymax = 0.0771;
%     sign_twin = -1;




% shotnum = [160127007 160127006;
%     160127009 160127008;
%     160127011 160127010;
%     160127014 160127012;
%     160127015 160127012;
%     160127017 160127016;
%     160127018 160127016;
%     160127020 160127019;
%     160127022 160127021;
%     160127023 160127021;
%     160127025 160127024;
%     160127026 160127024;
%     160127028 160127027;
%     160127029 160127027;
%     160127031 160127030;
%     160127032 160127030;
%     160127034 160127033;
%     160127035 160127033;
%     160127037 160127036;
%     160127038 160127036;
%     160127040 160127039;
%     160127041 160127039;
%     160127043 160127042;
%     160127044 160127042;
%     160127046 160127045;
%     160127047 160127045;
%     160127049 160127048;
%     160127050 160127048;
%     160127052 160127051;
%     160127053 160127051;
%     160127055 160127054;
%     160127056 160127054;
%     160127058 160127057;
%     160127059 160127057;
%     160127061 160127060;
%     160127062 160127060;
%     160127064 160127063;
%     160127065 160127063;
%     160127066 160127063;
%     160127068 160127067;
%     160127069 160127067;
% ];






%     shots = [151027024 151027023];


% for diagnosing case of equal side lengths:
%     shotnum_base = 150617012;
%     shotnum_def = 150617013;

% pulse for rev. of sci. instr. paper:
%     shots = [150915008 150915007];
%     shotnum_base = 150915007;
%     shotnum_def = 150915008;

% d =.40; %      % Reconstruction distance [m]
% xmin = 0.0058;  %0.006;
% xmax = 0.0188;  %0.019;
% ymin = 0.01685;  %0.015;
% ymax = 0.032;  %0.035;
% sign_twin = -1;


% nice pinch (shot w/ shadowgraphy effect):
%     shotnum_base = 151027023;
%     shotnum_def = 151027024;

%     shotnum_base = 151030026;
%     shotnum_def = 151030027;

%     shotnum_base = 151030030;
%     shotnum_def = 151030028;
% shots = [151027026 151027025;
%     151027040 151027041;
%     151027047 151027046;
%     151028011 151028012;
%     151028019 151028020;
%     151028027 151028025;
%     151028029 151028030;
%     151028031 151028032;
%     151028039 151028040;
%     151028060 151028061;
%     151028068 151028069;
%     151029012 151029014;
%     151029022 151029021;
%     151029028 151029029;
%     151029030 151029031];

%     shots =     [151030017 151030016;
%         151030024 151030022;
%         151030025 151030026;
%         151030027 151030026];
% [151027024 151027023;
%     151028034 151028035;
%     151028036 151028035;
%     151028053 151028052
%     151029015 151029014;
%     151029032 151029031;
%     151029033 151029034;
%     [151030017 151030016;
%     151030024 151030022;
%     151030025 151030026;
%     151030027 151030026];

% 151030000 + [   28 30;
%     28 26;
%     25 26;
%     25 22;
%     24 22;
%     23 22;
%     18 19;
%     18 16;
%     13 16];
% shots = 1;

% else
%     shots = NaN;
% end



%         shots = [160524000+[5:41]];
%         shots_tmp = shot_base_pair(shots,whichPair);
%         d =.5; %      % Reconstruction distance [m]
%         xmin = 0.0441;  %0.006;
%         xmax = 0.0569;  %0.019;
%         ymin = 0.0230;  %0.015;
%         ymax = 0.0425;  %0.035;
%         sign_twin = -1;
%         shots = shots_tmp;





%         shotnum = [160518000+[5:72]];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d =.50; %      % Reconstruction distance [m]
%         xmin = 0.0441;  %0.006;
%         xmax = 0.0569;  %0.019;
%         ymin = 0.0230;  %0.015;
%         ymax = 0.0425;  %0.035;
%         sign_twin = -1;
%         shots = shots_tmp;

%         shotnum = 160426000+[3:52];%103];
%         for i = 1:length(shotnum)-1
%             shots_tmp(i,:) = [shotnum(i+1) shotnum(i)];
%         end
% %         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(1:size(shots_tmp,1)) = 0.5;
%         xmin(1:size(shots_tmp,1)) = 0.0446;
%         xmax(1:size(shots_tmp,1)) = 0.0569;
%         ymin(1:size(shots_tmp,1)) = 0.0237;
%         ymax(1:size(shots_tmp,1)) = 0.0429;
%         sign_twin(1:size(shots_tmp,1)) = -1;
%         shots = shots_tmp;

%         shotnum = 160414000+[54:104];
%         for i = 1:length(shotnum)-1
%             shots_tmp(i,:) = [shotnum(i+1) shotnum(i)];
%         end
% %         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(1:size(shots_tmp,1)) = 1.2;
%         xmin(1:size(shots_tmp,1)) = 0.0911;
%         xmax(1:size(shots_tmp,1)) = 0.104;
%         ymin(1:size(shots_tmp,1)) = 0.0708;
%         ymax(1:size(shots_tmp,1)) = 0.0914;
%         sign_twin(1:size(shots_tmp,1)) = -1;
%         shots = shots_tmp;



%         shotnum = 160413000+[3:52];
%         for i = 1:length(shotnum)-1
%             shots_tmp(i,:) = [shotnum(i+1) shotnum(i)];
%         end
% %         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(1:size(shots_tmp,1)) = 1.1;
%         xmin(1:size(shots_tmp,1)) = 0.0843;
%         xmax(1:size(shots_tmp,1)) = 0.0973;
%         ymin(1:size(shots_tmp,1)) = 0.0667;
%         ymax(1:size(shots_tmp,1)) = 0.0869;
%         sign_twin(1:size(shots_tmp,1)) = -1;
%         shots = shots_tmp;

%         shotnum = 160323000+[7:50];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(1:size(shots_tmp,1)) = 1.2;
%         xmin(1:size(shots_tmp,1)) = 0.092;
%         xmax(1:size(shots_tmp,1)) = 0.105;
%         ymin(1:size(shots_tmp,1)) = 0.0737;
%         ymax(1:size(shots_tmp,1)) = 0.0943;
%         sign_twin(1:size(shots_tmp,1)) = -1;
%         shots = shots_tmp;
%
%         clear shotnum shots_tmp
%         shotnum = 160316000+[12:45];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 1.2;
%         xmin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.092;
%         xmax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.105;
%         ymin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0737;
%         ymax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0943;
%         sign_twin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = -1;
%         shots = [shots; shots_tmp];
%
%         clear shotnum shots_tmp
%         shotnum = 160224000+[7:55];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.85;
%         xmin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0671;
%         xmax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0799;
%         ymin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0538;
%         ymax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0740;
%         sign_twin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = -1;
%         shots = [shots; shots_tmp];
%
%         clear shotnum shots_tmp
%         shotnum = 160203000+[3:41];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.8128;
%         xmin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0639;
%         xmax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0769;
%         ymin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0508;
%         ymax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0711;
%         sign_twin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = -1;
%         shots = [shots; shots_tmp];
%
%         clear shotnum shots_tmp
%         shotnum = 160202000+[5:30];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.8128;
%         xmin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0639;
%         xmax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0769;
%         ymin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0508;
%         ymax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0711;
%         sign_twin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = -1;
%         shots = [shots; shots_tmp];
%
%         clear shotnum shots_tmp
%         shotnum = 160127000+[8:70];
%         shots_tmp = shot_base_pair(shotnum,whichPair);
%         d(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.88;
%         xmin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0707;
%         xmax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0838;
%         ymin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0570;
%         ymax(size(shots,1):size(shots,1)+size(shots_tmp,1)) = 0.0771;
%         sign_twin(size(shots,1):size(shots,1)+size(shots_tmp,1)) = -1;
%         shots = [shots; shots_tmp];