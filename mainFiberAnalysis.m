
clear;
clc;
close all;

%% User Input
pxSizeXY = 140; %in nm
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

idx = contains({file2Analyze.name},data2Use);
idx2Stack = find(idx);
nImStacks = length(idx2Stack);
%preallocate memory
allData = struct('filename',[],'pores2D',[],'pores3D',[],'polVolume',[],...
    'poreVolume',[],'totVolume',[],'ratioPol',[],'ratioPores',[]);
allData(nImStacks).filename = [];
h = waitbar(0);
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
        
        
    end
end

%% Fiber analysis
switch dim
    case '2D'
        [fiberProps2D,skel] = Calc.getFiberProps2D(IM);
    case '3D'
        disp('Analysis is not done yet')
end


%% simple plot
skel2Plot = skel(:,:,1);
branchPoints =  bwmorph(skel2Plot,'branchpoints');

SE = strel('square',3);
dilatedBranchPoints = imdilate(branchPoints,SE);

skelShift = double(skel2Plot);
skelShift(skel2Plot==0) = 0.5;
plotBranches = skelShift-dilatedBranchPoints;
plotBranches(plotBranches<0) = 0;
figure
imagesc(plotBranches)

