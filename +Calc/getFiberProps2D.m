function [fiber2D,allSkel] = getFiberProps2D(bw,data)
    fiber2D = struct();
    fiber2D.connectivity = [];
    fiber2D.numberOfNodes = [];
    fiber2D.numberOfBranches = [];
    fiber2D.thicknessStats = [];
    fiber2D.branchLength = [];
    
    allSkel = zeros(size(bw));
    for i = 1:size(bw,3)
        currentIm = bw(:,:,i);
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
        
        fiberImage = zeros(size(currentIm));
        figure(1)
        imagesc(currentIm);
        hold on
        
        thickness1 = zeros(length(branchLength),1);
        thickness2 = thickness1;
        for j = 1:length(branchLength)
            currentFiber = branchLength(j).PixelIdxList;
%           
            maxDist = max(distanceMap(currentFiber));
            
            [fiberInd, fiberCoord] = findFiberCenter(currentFiber,fiberImage);
            
            [idx,xVec,yVec] = getLinePixel(fiberCoord,size(fiberImage),maxDist);
            figure(1)
            hold on
            plot(xVec,yVec,'k','Linewidth',1.5);
            
            %get line profile
            lineData = double(data(idx));
            
            %get xAxis
            xAxis = sqrt((yVec-yVec(1)).^2 + (xVec-xVec(1)).^2);
            
%             figure(2)
%             plot(xAxis,lineData)
            guess.sig = double(maxDist/2);
            guess.mu =double( median(xAxis));
            guess.minMaxDomain = double([xAxis(1) xAxis(end)]);
            
           [FitPar,Fit,~]=SimpleFitting.gauss1D(lineData,double(xAxis),guess);
           
           
           thickness1(j) = distanceMap(fiberInd(round(length(currentFiber)/2)));
           thickness2(j) = FitPar(1);


        end
        fiber2D.numberOfBranches = [[fiber2D.numberOfBranches] ;length(branchLength)];
        fiber2D.connectivity = [[fiber2D.connectivity] ;connectivity];
        fiber2D.numberOfNodes = [[fiber2D.numberOfNodes] ;numberOfNodes];
        fiber2D.thicknessStats = [[fiber2D.thicknessStats]; thickness1,thickness2;];
        fiber2D.branchLength    = [[fiber2D.branchLength]; [branchLength.Area]'];
        allSkel(:,:,i) = skel;
    end

end


function [fiberInd,fiberCoord] = findFiberCenter(currentFiber,fiberImage)
    disp('test')
    
    fiberImage(currentFiber) = 1;
    
    
    endLine = bwmorph(fiberImage,'Endpoints');
    
    [startLine] = find(endLine);
    [xStartLine,yStartLine] = ind2sub(size(fiberImage),startLine);
    fiber = zeros(length(currentFiber),1);
    
    for i = 1:length(currentFiber)
        fiber(i) = [startLine(1)] ;
               
        %find which neighbor pixel is equal to one
        neighbor = findNeighbor([xStartLine(1) yStartLine(1)],size(fiberImage),1,8);
        
        idx = sub2ind(size(fiberImage),neighbor(:,1),neighbor(:,2));
        
        nextPxIdx = and(ismember(idx,currentFiber), ~ismember(idx,fiber));
        nextPx = idx(nextPxIdx);
        
        startLine = nextPx;
        [xStartLine,yStartLine] = ind2sub(size(fiberImage),startLine);
    end
    
    [fibRow,fibCol] = ind2sub(size(fiberImage),fiber);
    
    fiberInd = fiber;
    
    fiberCoord = [fibRow,fibCol];  

end

function [neighbor] = findNeighbor(idx,dim,r,neigh)
    %function to find neighboring pixel in a given radius of a
    %central pixel given by idx. The function also makes sure that
    %we do not have indices outside the images.
    
    switch neigh
        case 8
            iIdx = idx(1)-r:idx(1)+r;
            jIdx = idx(2)-r:idx(2)+r;

            iIdx = iIdx(iIdx>=1 & iIdx<=dim(1));
            jIdx = jIdx(jIdx>=1 & jIdx<=dim(2));

            neighbor = combvec(iIdx,jIdx)';
        case 4
            neighbor = [idx(1)-r,idx(2);
                        idx(1)+r, idx(2);
                        idx(1),idx(2);
                        idx(1),idx(2)-r;
                        idx(1),idx(2)+r;
                        ];
            
            idx2Del1 = neighbor(:,1)>=1 & neighbor(:,1)<=dim(1);
            idx2Del2 = neighbor(:,2)>=1 & neighbor(:,2)<=dim(2);
            
            idx2Del = and(idx2Del1,idx2Del2);
            
            neighbor(~idx2Del,:) = [];
            
        otherwise
            error('inconsistent number of neighbor chosen, please use 4 or 8')
    end
end



 function [idx,xVec,yVec] = getLinePixel(fiberCoord,dim,vecLength)
    idxCenter = round(length(fiberCoord)/2);
    centerPoint = fiberCoord(round(length(fiberCoord)/2),:);
    
    x = [fiberCoord(idxCenter+1,2); fiberCoord(idxCenter-1,2)];
    y = [fiberCoord(idxCenter+1,1); fiberCoord(idxCenter-1,1)];

    center = [mean(x),mean(y)];
    %get the slope of the line
    m = (y(2)-y(1))/(x(2)-x(1));
    
    % we want the perpendicular:
    mPerp = -1/m;
    %get the origin of the line
    b = centerPoint(1)-mPerp*centerPoint(2);
    %sort data point base on x axis
    [x,idx] = sort(x);
        
    if abs(mPerp) ~= inf
        xVec = centerPoint(2)-(ceil(vecLength/2+4)):centerPoint(2)+(ceil(vecLength/2+4));
        yVec = mPerp*xVec+b;
        
    else
        yVec = centerPoint(1)-(ceil(vecLength/2+4)):centerPoint(1)+(ceil(vecLength/2+4));
        xVec = ones(size(yVec))*x(1);

    end

    xVec = round(xVec);
    yVec = round(yVec);

    [xVec,yVec] = fixVec(xVec,yVec,dim);
    %get the index of the pixel in the image coordinates
    idx = sub2ind(dim,yVec,xVec);

    end

    function [xVec,yVec] = fixVec(xVec,yVec,dim)

        %delete value below 1
        idx1= xVec<1;
        idx2 = xVec>dim(2);
        idx3 = yVec<1;
        idx4 = yVec>dim(1);

        idx = idx1+idx2+idx3+idx4;
        idx(idx>1) = 1;
        idx = logical(idx);
        if ~all(idx==0)
            %delete value above dim
            xVec(idx) =[];
            yVec(idx) = [];
        end
    end
