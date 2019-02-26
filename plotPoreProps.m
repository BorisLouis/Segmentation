clear;
clc;
close all;
%% user input
nBins = 100;
dim   = '3D';
prop2Plot = 'connect';
fileExt = '.mat';
outputName = 'Figures';
%% Load 
[file2Analyze,currentFolderName,outDir] = Load.Folder(fileExt,outputName);
data2Plot = Load.getData2Plot(file2Analyze,dim,prop2Plot);
fields = fieldnames(data2Plot);
nFile  = length(fields);
%% GetCDF
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