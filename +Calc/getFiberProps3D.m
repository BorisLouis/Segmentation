function [fiber3D,skel] = getFiberProps3D(IM)
    
    fiber3D = struct();
    fiber3D.connectivity = [];
    fiber3D.numberOfNodes = [];
    fiber3D.numberOfBranches = [];
    fiber3D.thicknessStats = [];
    allSkel = zeros(size(IM));
   
    skel = bwskel(IM,'MinBranchLength',10);
    branchPoints =  bwmorph(skel,'branchpoints');

    %dilate branchpoints to provide clearcuts (not ideal for very small fibers)
    SE = strel('square',3);
    dilatedBranchPoints = imdilate(branchPoints,SE);

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

    fiber3D.connectivity = [[fiber3D.connectivity] ;connectivity];
    fiber3D.numberOfNodes = [[fiber3D.numberOfNodes] ;numberOfNodes];
    fiber3D.thicknessStats = [[fiber3D.thicknessStats]; thicknessStats];

end