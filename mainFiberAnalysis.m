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
outputName = 'Fiber-Results';
data2Use   = '_invAdapt';%_invAdapt %'_adapt'

%% Loading Data
pxArea = pxSize.XY*pxSize.XY; %in µm^2
pxVol  = pxSize.XY*pxSize.XY*pxSize.Z;%!!Image is rescaled before calculation

[path] = uigetdir();
mainFolderContent = dir(path);
mainFolderContent(~[mainFolderContent.isdir]) = [];
for i = 3:length(mainFolderContent)%avoid the . and ..
    if strcmp(mainFolderContent(i).name, 'SegmentedStacks')
        folder2Binary = [mainFolderContent(i).folder filesep mainFolderContent(i).name];
        folder2Data   = split([mainFolderContent(i).folder filesep mainFolderContent(i).name],'SegmentedStacks');
        folder2Data = folder2Data{1};
    else
        folder2Data = [mainFolderContent(i).folder filesep mainFolderContent(i).name];
        folder2Binary = [mainFolderContent(i).folder filesep mainFolderContent(i).name filesep 'SegmentedStacks'];
    end
    %get binary
    currentFolderContent = dir(folder2Binary);
    index2Images   = contains({currentFolderContent.name},fileExt);
    tmp = currentFolderContent(index2Images);
    idx = contains({tmp.name},data2Use);
    file2Analyze(i) = tmp(idx);
    %get Data
    dataFolderContent = dir(folder2Data);
    index2Images   = contains({dataFolderContent.name},fileExt);
    data2Analyze(i) = dataFolderContent(index2Images);
   

end
assert(~isempty(file2Analyze), sprintf('no %s found in the directory', fileExt));
assert(~isempty(data2Analyze), sprintf('no %s found in the directory', fileExt));
idx = or(strcmp('.', {mainFolderContent.name}),strcmp('..', {mainFolderContent.name}));
file2Analyze(idx) = [];        
data2Analyze(idx) = [];
%% Looping through the Data

h = waitbar(0);

nImStacks = length(file2Analyze);
%preallocate memory
allData = struct('filename',[],'fiber2D',[],'fiber3D',[]);
allData(nImStacks).filename = [];

%loop through stack
for i = 1:nImStacks
    
    hMessage = sprintf('Loading image stack number %d/%d',i,nImStacks);
    waitbar(i/nImStacks,h,hMessage);

    %Data loading
    path2Binary = strcat(file2Analyze(i).folder,filesep);
    
    tmpName = file2Analyze(i).name;
    %Check which data to be used
    if isempty(strfind(file2Analyze(i).name,data2Use))
    
    else
        
        p2file      = strcat(path2Binary,tmpName);
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
        binaryData = Load.Movie.tif.getframes(p2file,frames);
        
        % load Data
        path2Data = [data2Analyze(i).folder,filesep data2Analyze(i).name]; 
        fileInfo    = Load.Movie.tif.getinfo(path2Data);
        frames = 1:fileInfo.Frame_n;
        data = Load.Movie.tif.getframes(path2Data,frames);

        switch dim
            case '2D'
                
                %test correlation 
                [corrMat] = normxcorr2(data(:,:,1),data(:,:,1));
                [a,b] = Calc.radial_profile(corrMat,1);
                
                %test fft
                [dataFFT] = fftshift(fft2(data(:,:,1)));
                [a1,b1] = Calc.radial_profile(real(dataFFT),1);
                
                [fiberProps2D,skel] = Calc.getFiberProps2D(binaryData,data);
                
                %convert to micrometer
                fiberProps2D.branchLength = fiberProps2D.branchLength * pxSizeXY;
                fiberProps2D.thicknessStats = fiberProps2D.thicknessStats* pxSizeXY *2; %convert to diameter
                fiberProps2D.corrDecayX = a;
                fiberProps2D.corrDecayY = b;
                fiberProps2D.fftX = a1;
                fiberProps2D.fftY = b1;
                fiberProps3D = [];
            case '3D'
                error('Analysis is not done yet')
        end
        
        allData(i).filename = file2Analyze(i).name;
        allData(i).fiber2D = fiberProps2D;
        allData(i).fiber3D = fiberProps3D;
        
        
        filename = [path2Binary filesep 'skeleton.mat'];
        save(filename,'skel')

        filename = [path2Binary filesep 'skeleton.png'];
        f = figure; 
        imagesc(skel(:,:,1))
        axis image
        colormap('gray')
        saveas(gcf,filename,'png')
        
        close(f);

        
        
    end
end
%% save data
filename = [path filesep 'fiberProps.mat'];
save(filename,'fiberProps2D')


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
axis image
colormap('gray')

