% The aim of this code is to plot CCDF and histograme for the properties 
% previously calculated by mainPoreSizeCalc.m. It was coded to be
% relatively flexible so one can plot any property extracted

% HOW TO USE:
% In loading input the user provide information about the folder in which
% the data he wants to process is (typically PoreSize-Results if one is
% using our workflow.
% Then the user needs to provide info about the dimension of the analysis
% results he/she wants to plot (2D or 3D)and the pore property (volume,
% area, extRad, inRad, throats,...)

% The code will output a figure with 2 subplot, a CCDF and an Histogram and
% one with... (ADD SUSANA's CODE)

% Note: The order of the section allow the user to change the property
% to be plotted without having to load the data again by running the code
% section per section (shift+enter on the targeted section).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %
% Susana Rocha (https://github.com/SusanaRocha)                           %
%                                                                         %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %
%           Boris Louis:  https://borislouis.github.io/                   %
%           Susana Rocha: https://susanarocha.github.io/                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear;
clc;
close all;
%% Loading input
fileExt = '.mat';
outputName = 'Figures';
[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);
%% user input
nBins = 100;
dim   = '3D';
prop2Plot = 'connect';
%% Prepare Data to be plotted
data2Plot = Load.getData2Plot(file2Analyze,dim,prop2Plot);
fields = fieldnames(data2Plot);
nFile  = length(fields);
%% Get the CCDF
for i = 1 : nFile
    
    currField = fields{i};
    allData = [];
    
    for j = 1:size(data2Plot.(currField),1)
        currData  = data2Plot.(currField).(prop2Plot){j,:};
        [~,CCDF] = Plotting.getCDF(currData(:));
        data2Plot.(currField).CDF{j,1} = [CCDF.x CCDF.y];
        allData = [allData; currData(:)]; 
    end
         
    [~,CCDF] = Plotting.getCDF(allData);
    data2Plot.(currField).allCCDF{1,1} = [CCDF.x CCDF.y];
    data2Plot.(currField).allData{1,1} = allData;


end
%% plotting
%histogram

for i = 1 : nFile
    currField = fields{i};
    currData  = data2Plot.(currField);
    
    figure(i)
    subplot(1,2,1)
    histogram(currData.(prop2Plot){1,1})
    title(['Histogram of ' prop2Plot])
    axis square
    xlabel(prop2Plot);
    ylabel('Occurence')
    
    subplot(1,2,2)
    plot(currData.CDF{1,1}(:,1),currData.CDF{1,1}(:,2))
    ax = gca;
    title(['CCDF' prop2Plot])
    ax.XScale = 'log';
    ax.YScale = 'log';
    axis square
    xlabel([prop2Plot ' Log scale']);
    ylabel('[1-CDF] Probability [0 1]')
    
end