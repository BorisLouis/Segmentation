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
        disp('stop')
        
    end
end


%% test skeletonize
tmp = IM;
% THIS IS BAD METHOD
%skel = bwskel(tmp,'MinBranchLength',2);
%skel = bwmorph(tmp,'skel',inf);.

skel = bwskel(tmp,'MinBranchLength',10);

figure
subplot(1,2,1)
imagesc(labeloverlay(double(tmp),double(skel),'Transparency',0))
subplot(1,2,2)
imagesc(labeloverlay(double(tmp),double(skel),'Transparency',0))

%% Calculate skeleton properties 


branchPoints =  bwmorph(skel,'branchpoints');
%dilate branchpoints
SE = strel('disk',2);
dilatedBranchPoints =imdilate(branchPoints,SE);

figure
imagesc(skel-dilatedBranchPoints)

%D = bwdistgeodesic(skel2,find(B),'quasi');
%for plotting:
D2 = graydist(skel,find(branchPoints));

branches = skel-dilatedBranchPoints;
branches(branches<0) = 0;

labelSkel = bwlabel(branches);
figure
imagesc(labelSkel)

branchLength = regionprops(labelSkel,'Area');

labelledNode = bwlabel(branchPoints);
numberOfNodes = max(labelledNode(:));

% connectivity 
%==> same pore connectivity

%thickness

% distMap = bwdist(~IM);
% ws = watershed(-distMap);
% ws(~IM) = 0;
% 
% figure
% imagesc(labeloverlay(double(ws),double(branches),'Transparency',0))


%% test from internet
inverse = imcomplement(IM);
ED = edge(IM,'Sobel');


BP = bwmorph(skel, 'branchpoints');
EP = bwmorph(skel, 'endpoints');
[y,x] = find(EP);
BP_L = find(BP);
Dmask = false(size(skel));

%D = bwdistgeodesic(skel);
% % tic
% for k = 1:length(x)
%     D = bwdistgeodesic(skel,x(k),y(k));
%     distanceToBranchPt = min(D(BP_L));
%     Dmask(D < distanceToBranchPt) = true;
% end
% toc

% finalSkel = skel - Dmask;
figure(1)
muxing = ED | inverse;
imshowpair(muxing,skel,'blend')

%% thickness
finalIm = muxing+skel;

distance = bwdist(finalIm);
figure
imagesc(distance);


skel = bwskel(distance,'MinBranchLength',10);


%%



test = double(IM);
test(test==0) = -inf;

figure
imagesc(test)

test = test-branches;

figure
imagesc(~test)
