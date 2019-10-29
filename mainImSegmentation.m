% The aim of this code is to analyze a folder containing Z-stack with    
% structure of interest. In our case we used it for characterizing       
% porous material. 
%
% The program segments the stacks (e.g. binarization pore-material) and
% save the information for further analysis
% HOW TO USE:
% ==> User input can be modified and are described in the code section
% below
% ==> The code will open the file selecter from window, A test file is 
% provided (testFile.tif) so you can give it a try.
% ==> output binary images are saved in the folder where the file was, in a
% new folder called "SegmentedStacks
% ==> plot of the raw data with the resulting segmentation is shown, the
% user can switch from frame to frame by clicking
%
% !!!!!!!!!!!!!!! The output data is a reverse binary image (what was
% considered dark on the image (e.g. pores) will be bright (1) on the 
% binary data while what was considered as material (bright on the image)
% will be dark (0) on the binary data. This was done because then
% mainPoreSizeCalc can directly use the data to calculate the pore 
% properties
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %
%                                                                         %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %
%           Boris Louis: https://borislouis.github.io/                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear 
close all
clc
%% User Input
threshold = 0.2; %sensitivity for adaptive threshold
connectivity = 216; %3D connectivity for binarization (only for adaptive threshold)
diskDim = 4; %disk dimension to clean segmentation artefact (the bigger the more get cleaned).
S = 2; % size of gauss filter in pixel
pixZ  = 2;% size of pixel in z vs x/y (ratio)
fileExt = '.tif'; % only work with TIF
toAnalyze = 'file'; %if folder is chosen the code will search for all the tif
%inside the folder and segment them with the input provided.
%if 'file' is input only one file will be segmented.
memEff = true; % memory efficiency. For large files the code struggle doing
%the segmentation on standard computer, thus we encoded a possibility to
%process the stack in several steps of dIDX x dIDX (which is a standard format
%for microscopy camera. The steps results are then put together, this can 
%created some small artefact at the interface of these pieces but was mostly
%not an issue in our case.
dIDX = 512; % to split stack in 512*512*Z chunks
%% Loading Data
outputName  = 'SegmentedStacks';%folder name where output stack will be saved
switch toAnalyze
    case 'folder'
        %Load folder, and create a folder for data output.
        [file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);
        assert(~isempty(file2Analyze), sprintf('no %s found in the directory', fileExt));
    case 'file'
        [fileName, path] = uigetfile('*.tif', 'Select a file to segment');
        file2Analyze.name = fileName;
        file2Analyze.folder = path;
        [~,~,ext] = fileparts(fileName);
        outDir = [path filesep outputName];
        status = mkdir(outDir); 
    otherwise
        error('Unknown type to analyze');
end

%% Processing
nFiles = size(file2Analyze,1);

for i = 1 : nFiles

    disp(['Loading stack --------------' file2Analyze(i).name])
    %getting file info
    path2Stacks = strcat(file2Analyze(i).folder,filesep);
    tmpName = file2Analyze(i).name;
    p2file      = strcat(path2Stacks,tmpName);
    fileInfo    = Load.Movie.tif.getinfo(p2file);
    
    %remove warnings for loading tif
    warning('off','all');
    
    %Check number of Frame
    tNframes = fileInfo.Frame_n;
    frames2load = 1:tNframes;
    %loading occurs here
    IM     = Load.Movie.tif.getframes(p2file, frames2load); %Loading on of the frame
    warning('on','all');
    disp('DONE with loading --------------')

    %%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%
    disp('Now doing 3D gauss filtering this can take about 3 minutes')
   
    sigma = [S,S,S/pixZ];
    IMs = imgaussfilt3(IM, sigma);%3D gaussian filtering
    disp('DONE with filtering ------------')

    %%%%%%%%%%%%%%% Segmenting %%%%%%%%%%%%%%%
    disp('Now doing segmentation this can take a few minutes ~10')
    imSize = size(IMs);
    
    if imSize(1)<dIDX || imSize(2)<dIDX
        memEff = false;
        warning('cannot cut into chunks as data is smaller that the requested chunk size');
    end
    
    if memEff
        %pre allocate memory for storing binary images
        BWglobal = false(imSize);
        BWadapt  = false(imSize);
        % get list of indices
        mod1 = mod(imSize(1),dIDX);
        mod2 = mod(imSize(2),dIDX);
        assert(or(mod1==0,mod1>50), 'Your chunck size is generating problems on 1st index')
        assert(or(mod2==0,mod2>50), 'Your chunck size is generating problems on 2nd index')

        xi = 1:dIDX:imSize(1);
        xf = dIDX:dIDX:imSize(1);
        if xf(end)~= imSize(1)
            xf = [xf, imSize(1)];
        end
        assert(length(xi)==length(xf), 'Problems! indexing over the volume')
        % get the position of the different stack
        IDX = zeros(length(xi)*length(xf),4);
        [XX, YY] = meshgrid(xi,xf);
        IDX(:,1) = XX(:);
        IDX(:,4) = YY(:);
        [XX, YY] = meshgrid(xf,xi);
        IDX(:,2) = XX(:);
        IDX(:,3) = YY(:);
        clear XX YY xi xf
        
        %get the number of chunks
        lastIDX = size(IDX,1);
        lastStr = num2str(lastIDX);
     
        h = waitbar(0, sprintf('Starting segmentation of stack %d/%d...',i,nFiles));
        for j = 1:lastIDX
            
            iStr = num2str(j);
            xi = IDX(j,1);
            xf = IDX(j,2);
            yi = IDX(j,3);
            yf = IDX(j,4);
            %extract current chunck from the data
            tmp = IMs(xi:xf,yi:yf,:);
            %segmentation occurs here
            [gBW,aBW] = imSegmentation.segmentStack(tmp,'threshold',threshold,...
                'connectivity',connectivity,'diskDim',diskDim);
            %store the segmented chunk at the right place in the binary
            %stack
            BWadapt (xi:xf,yi:yf,:) = aBW;%adaptive threshold
            disp(['Done adaptive step ' iStr '/' lastStr ])
            BWglobal(xi:xf,yi:yf,:) = gBW;%global threshold
            disp(['Done global step ' iStr '/' lastStr ])
            
            waitbar(j/lastIDX,h,sprintf('Segmentation of stack %d/%d... %d percent achieved',...
                i,nFiles, round(j/lastIDX*100)));
        end
        close(h);
    else
        %segmentation occurs here if single file was analyzed
        [BWglobal,BWadapt] = imSegmentation.segmentStack(IMs,'threshold',threshold,...
                'connectivity',connectivity,'diskDim',diskDim);
            dIDX = imSize(1);
    end

    %%%%%%%%%%%%%%% Data saving %%%%%%%%%%%%%%%
    
    % now we save segmented images
    disp('Storing global results')
    tifName = [outDir filesep 'Seg_global_' tmpName];
    dataStorage.BinaryTiff(tifName,BWglobal);

    disp('Storing adaptive results')
    tifName = [outDir filesep 'Seg_adapt_' tmpName];
    dataStorage.BinaryTiff(tifName,BWadapt);
    
    %Save useful information about how the binarization was done
    infoFileName = [outDir filesep 'info.txt'];
    fid = fopen(infoFileName,'wt');
    fprintf(fid,'This file contains information intended for the user of segmentStack.m\n');
    fprintf(fid,' In such a way that the user knows what variable value were used.\n\n');
    fprintf(fid,'Adaptive Threshold sensitivity: %0.1f\n',threshold);
    fprintf(fid,'Number of frame analyzed: %d/%d\n',tNframes);
    fprintf(fid,'3D Gaussian filtering: S = %d; pixZ  = %d; sigma = [S,S,S/pixZ]\n',S,pixZ);
    fprintf(fid,'Three-dimensional connectivity - BWareaopen: %d\n',connectivity);
    fprintf(fid,'strel disk dimension for imopen: %d\n',diskDim);
    fprintf(fid,'Image partition in chunck: %d x %d x %d \n',dIDX,dIDX,tNframes);
    fclose(fid);
   
end
%output a message to tell that everything went ok
h = msgbox('The Data were succesfully saved !', 'Success');
%% Plotting the segmentation to check 
idx = 1;
path = file2Analyze.folder;
frameSkip = 10;%number of frame to skip between two clicks when going through the stack
imSegmentation.check(path,idx,frameSkip);

