%% Evaluate segmentation
% The aim of this code is to evaluate the performance of automated 
% segmentation method by comparing it with reference segmentation 
% (typically hand segmentation), it calculates precision, recall and
% bfscore

% The program segments load segmented stack and compare them 
% HOW TO USE:
% ==> Path of the reference segmentation needs to be provided in the user
% input part of the code
% ==> Press run
% ==> The code will open the file selecter from window, you need to select
% a segmented stack.
% ==> Gives out precision recall and bfscore in the command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Website : Boris Louis: https://borislouis.github.io/                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear;
close all;
clc;
%% User Input

path2HandSegmentation = 'D:\Documents\Unif\PhD\Papers\07 - PIC structural characterization\Revision paper\Sensitivity test\HandSegVolume\0_25\Hand\Manually-SegmentedStackStack0_25mgmL.tif';
%%

[file, path] = uigetfile('*.tif');
path2Seg = [path filesep file];

%load Ground truth segmentation
fileInfo    = Load.Movie.tif.getinfo(path2HandSegmentation);
warning('on','all')
tNframes    = fileInfo.Frame_n;
frames      = 1:fileInfo.Frame_n;
groundTruth = Load.Movie.tif.getframes(path2HandSegmentation,frames);

%load software segmented stack
fileInfo    = Load.Movie.tif.getinfo(path2Seg);
warning('on','all')
tNframes    = fileInfo.Frame_n;
frames      = 1:fileInfo.Frame_n;
segData     = Load.Movie.tif.getframes(path2Seg,frames);

%% Calculate evalutation metrics

[score, precision, recall] = bfscore(logical(segData),~logical(groundTruth));

disp(['The segmentation precision is ' num2str(precision)]);
disp(['The segmentation recall is ' num2str(recall)]);
disp(['The segmentation score is ' num2str(score)]);


y = [precision; recall; score];
b = bar(y);
set(gca,'xTickLabel',{'precision','recall','bfscore'})


