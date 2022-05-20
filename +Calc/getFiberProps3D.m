function [fiber2D,skel] = getFiberProps3D(IM)
    error('Work in progress')
    fiber2D = struct();
    fiber2D.connectivity = [];
    fiber2D.numberOfNodes = [];
    fiber2D.numberOfBranches = [];
    fiber2D.thicknessStats = [];
    allSkel = zeros(size(IM));
   
    skel = bwskel(currentIm,'MinBranchLength',10);
    branchPoints =  bwmorph(skel,'branchpoints');

    %dilate branchpoints to provide clearcuts (not ideal for very small fibers)
    SE = strel('square',3);
    dilatedBranchPoints = imdilate(branchPoints,SE);

%         skelShift = double(skel);
%         skelShift(skel==0) = 0.5;
%         plotBranches = skelShift-dilatedBranchPoints;
%         plotBranches(plotBranches<0) = 0;
%         figure
%         imagesc(plotBranches)
    % extract the branches by subtracting the branchpoints
    branches = skel-dilatedBranchPoints;
    branches(branches<0) = 0;

    labelSkel = bwlabel(branches);
%         figure
%         imagesc(labelSkel)

    branchLength = regionprops(labelSkel,'Area','PixelIdxList');
    idx = [branchLength.Area]<3;
    branchLength(idx) = [];

    labelledNode = bwlabel(dilatedBranchPoints);
    numberOfNodes = max(labelledNode(:));

    % connectivity 
    connectivity = zeros(numberOfNodes,1);

    SE = strel('square',3);
    for j = 1:numberOfNodes
        currentNode = labelledNode==j;

        dilatedNode = imdilate(currentNode,SE);

        overlap = labelSkel;
        overlap(~dilatedNode) = 0;

        %find overlap between fiber and current Node
        connectivity(j) = length(unique(overlap))-1;

    end

    %thickness
    inverse = imcomplement(currentIm);
    distanceMap = bwdist(inverse);

    %iterate through branches to get thickness
    thicknessStats = table(cell(length(branchLength),1),zeros(length(branchLength),1),...
        zeros(length(branchLength),1),'VariableNames',{'Data','Median','Std'});

    for j = 1:length(branchLength)
        currentFiber = branchLength(j).PixelIdxList;

        thicknessStats.Data{j} = distanceMap(currentFiber);
        thicknessStats(j,:).Median = median(thicknessStats.Data{j});
        thicknessStats(j,:).Std = std(thicknessStats.Data{j});


    end

    fiber2D.connectivity = [[fiber2D.connectivity] ;connectivity];
    fiber2D.numberOfNodes = [[fiber2D.numberOfNodes] ;numberOfNodes];
    fiber2D.thicknessStats = [[fiber2D.thicknessStats]; thicknessStats];

end