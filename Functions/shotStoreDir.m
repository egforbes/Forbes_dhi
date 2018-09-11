function [directory] = shotStoreDir(shotnum) 
% shotnum = 101021010;


homeBase = pwd;
SN = num2str(shotnum);

cd H:\Reconstructions
% dos(pushd \\zappc_1\Holography)


dateVect = [2000+str2num(SN(1:2)), ...
    str2num(SN(3:4)), str2num(SN(5:6)),0,0,0];
dateS = datestr(dateVect);
year = dateS(8:end); month = dateS(4:6); day = dateS(1:2);
DM = [day,month];
currentDirY = dir;
dirNamesY = {currentDirY.name};
if ~strcmp(dirNamesY,year)
    mkdir(year);
end

cd(year)
currentDirDM = dir;

dirNamesDM = {currentDirDM.name};

if ~strcmp(dirNamesDM,DM)
    mkdir(DM)
end

cd(DM)

warning off;

mkdir('phase_diff')
% mkdir('fig1')
% mkdir('fig2')
% mkdir('fig3')
% mkdir('fig4')
% mkdir('fig5')
% mkdir('fig6')
% mkdir('fig7')


warning on;

cd(homeBase)

directory = ['H:\Reconstructions','\',year,'\',DM];