clear;
clc;
close all;
%% User input
path = 'F:\TIF files';
folder2Plot = {'Cold 1', 'Hot 1', 'Cold 2'};


%% Data Loading
mainFolder = dir(path);
data2Plot = struct();
dataName = 'fiberProps.mat';
for i = 1:length(folder2Plot)
    
    currentFolder2Plot = folder2Plot{i};
    idx = contains({mainFolder.name},currentFolder2Plot);
    
    data2Load = [path filesep mainFolder(idx).name filesep dataName];
    
    tmpData = load(data2Load);
    
    varName = strrep(currentFolder2Plot,' ','');
    data2Plot.(varName) = tmpData.allData;

    
end
%% Pre-processing
avgData2Plot = struct();

for i = 1:length(folder2Plot)
    currentFolder2Plot = folder2Plot{i};
    varName = strrep(currentFolder2Plot,' ','');
    
    currentData = data2Plot.(varName);
    allData = [];
    for j = 1:length(currentData)
        tmpData = currentData(j).fiber2D;
        fields = fieldnames(tmpData);
        
        for k = 1:length(fields)
            currentField = fields{k};
            if strcmp(currentField,'numberOfNodes') || strcmp(currentField,'numberOfBranches')...
                    || strcmp(currentField,'thicknessStats') || strcmp(currentField,'branchLength')...
                    || strcmp(currentField,'straightness') || strcmp(currentField,'ratio')...
                    || strcmp(currentField,'porosity')
                if j == 1
                    allData.(currentField) = tmpData.(currentField);
                else
                    allData.(currentField) = [allData.(currentField); tmpData.(currentField)];
                end
            end
            
            if strcmp(currentField, 'IntProb') || strcmp(currentField, 'IntFiber')
                x = 1:255;
                if j == 1
                    allData.Int = x;
                    tmpProb = zeros(size(x));
                    intProb = tmpData.(currentField);
                    int     = tmpData.Int;
                    tmpProb(ismember(int,x)) = intProb;
                    allData.(currentField) = tmpProb;
            
                else
                    tmpProb = zeros(size(x));
                    intProb = tmpData.(currentField);
                    int     = tmpData.Int;
                    tmpProb(ismember(int,x)) = intProb;
                    allData.(currentField) = allData.(currentField) + tmpProb;
                    
                end
            
            end

        end            
    end
    %normalize probability of observing a specific intensity.
    allData.IntProb = allData.IntProb/length(currentData);
    
    avgData2Plot.(varName) = allData;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING OCCUR HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fields = fieldnames(avgData2Plot);
%%  Intensity Plot
figure
hold on
for i = 1:length(fields)
   currentDataX = avgData2Plot.(fields{i}).Int;
   currentDataY =  avgData2Plot.(fields{i}).IntProb;
   
   plot(currentDataX,currentDataY)
    
end
axis square
box on
legend(fields)
xlabel(['Intensity(a.u.)'])
ylabel(['Probability'])

title('Intensity distribution')

%%  Intensity Fiber Plot
fields = fieldnames(avgData2Plot);

figure
hold on
for i = 1:length(fields)
   currentDataX = avgData2Plot.(fields{i}).Int;
   currentDataY =  avgData2Plot.(fields{i}).IntFiber;
   
   plot(currentDataX,currentDataY)
    
end
axis square
box on
legend(fields)
xlabel(['Intensity(a.u.)'])
ylabel(['Probability'])
title('Segmented Intensity distribution')
%% Number of Nodes
pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).numberOfNodes;
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).numberOfNodes)) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Number of Nodes')
title('Number of Nodes')

%% Number of branch
pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).numberOfBranches;
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).numberOfBranches)) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Number of Branches')
title('Number of Branches')
%% Thickness from distance map
pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).thicknessStats(:,1);
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).thicknessStats(:,1))) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Thickness - distmap')
title('Thickness - distmap')

%for Johannes to try log scale
%set(gca,'YScale','log')


   
figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of Thickness');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});


%% Thickness from line profile
pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).thicknessStats(:,2);
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).thicknessStats(:,2))) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Thickness - line prof.')
title('Thickness - line prof.')

figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of Thickness from line Prof');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});


%% branch length
 pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).branchLength;
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).branchLength)) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Branch length(nm)')
title('Branch length(nm)')

figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of branch length');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});
%% straightness
pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).straightness;
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).straightness)) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Straightness')
title('Straightness')

figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of straightness');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});

%% ratio from distMap
 pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).ratio(:,1);
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).ratio(:,1))) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Length to thickness ratio - distmap')
title('Length to thickness ratio - distmap')

figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of length to thickness ratio');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});
%% ratio from line profile
 pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).ratio(:,2);
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).ratio(:,2))) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Length to thickness ratio - line prof.')
title('Length to thickness ratio - line prof.')

figure 
Plotting.distributionPlot(data,'histOpt',0, 'showMM',0);
pp=findobj(gca, 'type', 'patch');
set(pp(1), 'FaceColor', [0 0 1]);  
set(pp(2), 'FaceColor', [1 0 0]);
set(pp(3), 'FaceColor', [0 0 0.7]);
hold on
title('Violin dist of length to thickness ratio line prof');
axis square
box on
xlabel('Sample')
%ylabel([prop2Plot,' (µm)'])
set(gcf, 'color', 'w')
set(gca, 'XTick', 1, 'XTickLabel', {'Cold','Hot','Cold'});
%% Porosity
 pts ={};
catIdx = [];
for i = 1:length(fields)
    data{i} =  avgData2Plot.(fields{i}).porosity;
    catIdx(end+1:end+length(avgData2Plot.(fields{i}).porosity)) = i;
    
end

figure
Plotting.PlotSpread.plotSpread(data,'categoryIdx',catIdx,...
'categoryMarkers',{'o','o','o'},'categoryColors',{'b','r','b'},'showMM',1)
axis square
box on
xlabel('Cold Hot Cold')
ylabel('Porosity')
title('Porosity')


