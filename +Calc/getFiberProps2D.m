function [fiber2D,allSkel] = getFiberProps2D(IM)
    fiber2D = struct();
    fiber2D.connectivity = [];
    fiber2D.numberOfNodes = [];
    fiber2D.numberOfBranches = [];
    fiber2D.thicknessStats = [];
    fiber2D.branchLength = [];
    
    allSkel = zeros(size(IM));
    for i = 1:size(IM,3)
        currentIm = IM(:,:,i);
        skel = bwskel(currentIm,'MinBranchLength',10);
        branchPoints =  bwmorph(skel,'branchpoints');

        %dilate branchpoints to provide clearcuts (not ideal for very small fibers)
        SE = strel('square',3);
        dilatedBranchPoints = imdilate(branchPoints,SE);

        % extract the branches by subtracting the branchpoints
        branches = skel-dilatedBranchPoints;
        branches(branches<0) = 0;

        labelSkel = bwlabel(branches);

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
        fiber2D.numberOfBranches = [[fiber2D.numberOfBranches] ;length(branchLength)];
        fiber2D.connectivity = [[fiber2D.connectivity] ;connectivity];
        fiber2D.numberOfNodes = [[fiber2D.numberOfNodes] ;numberOfNodes];
        fiber2D.thicknessStats = [[fiber2D.thicknessStats]; thicknessStats];
        fiber2D.branchLength    = [[fiber2D.branchLength]; [branchLength.Area]'];
        allSkel(:,:,i) = skel;
    end

end
