% The aim of this code is to calculate the pore size and other properties 
% of pores based on a pre-segmented image.The code assumes that the 
% segmentation is performed in a way that the pore are bright (=1) and the
% strtucture of interest is black which is the way our image segmentation
% works.

% HOW TO USE:
% ==> User input can be modified and are described in the code section
% below
% ==> The code will open the folder selecter from window, the user should
% select a folder where segmented data can be found (e.g. SegmentedStacks)
% ==> output structure is saved in a new folder within SegmentedSTacks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %
%                                                                         %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %
%           Boris Louis: https://borislouis.github.io/                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear;
clc;
close all;

%% User Input
pxSizeXY = 180; %in nm
pxSizeZ = 404;
pxSize.XY = pxSizeXY*10^-3;
pxSize.Z  = pxSizeZ*10^-3;
dim = '2D';%dimension to perform analysis: 'bubble', '2D', '3D' or 'both'
fileExt = '.tif';
outputName = 'PoreSize-Results';
data2Use   = 'adapt';%'global' or %'adapt'

%% Loading Data
% conversion from pixel to area
pxArea = pxSize.XY*pxSize.XY; %in µm^2
pxVol  = pxSize.XY*pxSize.XY*pxSize.Z;%!!Image is rescaled before calculation
% load folder containing binary file
[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);

%% Looping through the Data

h = waitbar(0);

idx = contains({file2Analyze.name},data2Use);
idx2Stack = find(idx);
nImStacks = length(idx2Stack);
%preallocate memory
allData = struct('filename',[],'pores2D',[],'pores3D',[],'polVolume',[],...
    'poreVolume',[],'totVolume',[],'ratioPol',[],'ratioPores',[]);
allData(nImStacks).filename = [];

%loop through stack
for i = 1:length(idx2Stack)
    idx = idx2Stack(i);
    hMessage = sprintf('Loading image stack number %d/%d',i,nImStacks);
    waitbar(i/nImStacks,h,hMessage);

    %Data loading
    path2Stacks = strcat(file2Analyze(idx).folder,filesep);
    tmpName = file2Analyze(idx).name;
    %Check which data to be used
    if isempty(strfind(file2Analyze(i).name,data2Use))
    
    else
        
        p2file      = strcat(path2Stacks,tmpName);
        warning('off','all')
        fileInfo    = Load.Movie.tif.getinfo(p2file);

        warning('on','all')
        tNframes = fileInfo.Frame_n;

        % init data that contains all infor for a single tif file
        tifStackData = [];

        %loop through the frames of the current stack
        nIM = tNframes;
        % waitbar(j/nImStacks,h,hMessage);
        frames = 1:fileInfo.Frame_n;
        IM = Load.Movie.tif.getframes(p2file,frames);

        [bubbles,pores2D,pores3D] = Calc.getPoreProps(IM,dim,pxSize);
        if ~isempty(pores2D)
            pores2D.area = pores2D.area*pxArea;
            pores2D.inRad = pores2D.inRad*pxSizeXY*1e-3;%in um
            pores2D.extRad = pores2D.extRad*pxSizeXY*1e-3;
            pores2D.throats = pores2D.throats*pxSizeXY*1e-3;
            pores2D.diameter = pores2D.diameter*pxSizeXY*1e-3;
            pores2D.ctr_ext = pores2D.ctr_ext*pxSizeXY*1e-3;
        end
        if ~isempty(pores3D)
            pores3D.vol = pores3D.vol*pxVol;
            pores3D.extRad = pores3D.extRad;
            pores3D.inRad = pores3D.inRad;
            pores3D.throats = pores3D.throats;
            pores3D.ctr_ext = pores3D.ctr_ext.*[pxSize.XY,pxSize.XY,pxSize.Z];
        end
%         bubbles.rad = bubbles.rad*pxSizeXY*1e-3;
%         bubbles.coord = bubbles.coord*pxSizeXY*1e-3;
        totVol = numel(IM);
        polVolume  = sum(sum(sum(~IM)));%Assume pore are brigth
        disp('Storing Data')

        allData(i).filename = file2Analyze(i).name;
        allData(i).pores2D = pores2D;
        allData(i).pores3D = pores3D;
        allData(i).polVolume(1) = polVolume;%stillPX
        allData(i).poreVolume(1)  = totVol-polVolume;
        allData(i).totVolume(1)    = totVol;%stillPX
        allData(i).ratioPol(1) = (polVolume)/totVol;
        allData(i).ratioPores(1) = (totVol-polVolume)/totVol;
        allData(i).bubbles = bubbles;
        allData(i).pxSizeXY = pxSizeXY;
        allData(i).pxSizeZ  = pxSizeZ;
        
    end
        disp('---------------------NEXT TIF ----------')
end
close(h);

%% saving data
infoFileName = [outDir filesep 'info.txt'];
fid = fopen(infoFileName,'wt');
fprintf(fid,'This file contains information intended for the user of poreSizeCalc\n');
fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
fprintf(fid,'Pixel size used: %d',pxSizeXY);
fclose(fid);

save([ outDir filesep data2Use 'PoreProps'],'allData','-v7.3');
h = msgbox('The Data were succesfully saved !', 'Success');
