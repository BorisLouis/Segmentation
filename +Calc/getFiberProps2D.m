function [fiber2D,allSkel] = getFiberProps2D(bw,data)
    fiber2D = struct();
    fiber2D.connectivity = [];
    fiber2D.numberOfNodes = [];
    fiber2D.numberOfBranches = [];
    fiber2D.thicknessStats = [];
    fiber2D.branchLength = [];
    fiber2D.straightness = [];
    fiber2D.ratio = [];
    allSkel = zeros(size(bw));
    
    data = imgaussfilt(data);
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
        imagesc(imfuse(data,skel));
        axis image
        hold on
        
        thickness1 = zeros(length(branchLength),1);
        thickness2 = thickness1;
        fiberStraightness = thickness1;
       
        fiberRatio1 = thickness1;
        fiberRatio2 = thickness2;
        
        for j = 1:length(branchLength)
            currentFiber = branchLength(j).PixelIdxList;
            [fiberInd, fiberCoord] = findFiberCenter(currentFiber,fiberImage);
            
            if ~isempty(fiberInd)
                maxDist = distanceMap(fiberInd(round(length(currentFiber)/2)));
             %get fiber straightness   
                nodeDistance = sqrt((fiberCoord(1,1)-fiberCoord(end,1)).^2+(fiberCoord(1,2)-fiberCoord(end,2)).^2);
                pathIntegral = 0;
                for k = 1:length(fiberCoord)-1

                    pathIntegral = pathIntegral + sqrt((fiberCoord(k,1)-fiberCoord(k+1,1)).^2+(fiberCoord(k,2)-fiberCoord(k+1,2)).^2);

                end

                fiberStraightness(j) = nodeDistance/pathIntegral;       
                
                %fiber thickness
                [idx,xVec,yVec] = getLinePixel(fiberCoord,size(fiberImage),maxDist);
                
                bwLineData = double(bw(idx));
                testLine = diff(bwLineData);
                
                if sum(abs(testLine))>=2
                    figure(1)
                    hold on
                    plot(xVec,yVec,'w','Linewidth',1.5);
                    %get line profile
                    lineData = double(data(idx));

                    %get xAxis
                    xAxis = sqrt((yVec-yVec(1)).^2 + (xVec-xVec(1)).^2);
                    
                    guess.sig = double(maxDist);
                    guess.mu =double( median(xAxis));
                    guess.minMaxDomain = double([xAxis(1) xAxis(end)]);

                   [FitPar,Fit,~]=SimpleFitting.gauss1D(lineData,double(xAxis),guess);

                    thickness1(j) = distanceMap(fiberInd(round(length(currentFiber)/2)));
                    thickness2(j) = FitPar(1);
                 

                    fiberRatio1(j) = pathIntegral/thickness1(j);
                    fiberRatio2(j) = pathIntegral/thickness2(j);
               
                else
                    thickness1(j) = NaN;
                    thickness2(j) = NaN;
                    fiberRatio1(j) = NaN;
                    fiberRatio2(j) = NaN;
               
                end
            else
               thickness1(j) = NaN;
               thickness2(j) = NaN;
               fiberRatio1(j) = NaN;
               fiberRatio2(j) = NaN;
               
            end
            


        end
        fiber2D.numberOfBranches = [[fiber2D.numberOfBranches] ;length(branchLength)];
        fiber2D.connectivity = [[fiber2D.connectivity] ;connectivity];
        fiber2D.numberOfNodes = [[fiber2D.numberOfNodes] ;numberOfNodes];
        fiber2D.thicknessStats = [[fiber2D.thicknessStats]; thickness1,thickness2;];
        fiber2D.branchLength    = [[fiber2D.branchLength]; [branchLength.Area]'];
        fiber2D.straightness = [[fiber2D.straightness]; [fiberStraightness]];
        fiber2D.ratio = [[fiber2D.ratio]; fiberRatio1,fiberRatio2;];
        
        
        allSkel(:,:,i) = skel;
    end

end


function [fiberInd,fiberCoord] = findFiberCenter(currentFiber,fiberImage)
       
    fiberImage(currentFiber) = 1;
    
    
    endLine = bwmorph(fiberImage,'Endpoints');
    
    [startLine] = find(endLine);
    [xStartLine,yStartLine] = ind2sub(size(fiberImage),startLine);
    fiber = zeros(length(currentFiber),1);
    if ~isempty(startLine)
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
    else
        fiberInd = [];
        fiberCoord = [];
    end

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
    shift = 2;
    
    id = [idxCenter - shift idxCenter + shift];
    id(id<1) = 1;
    id(id>size(fiberCoord,1)) = size(fiberCoord,1);
    
    x = [fiberCoord(id(2),2); fiberCoord(id(1),2)];
    y = [fiberCoord(id(2),1); fiberCoord(id(1),1)];

    center = [mean(x),mean(y)];
    %get the slope of the line
    m = (y(2)-y(1))/(x(2)-x(1));
    
    % we want the perpendicular:
    mPerp = -1/m;
    %get the origin of the line
    b = centerPoint(1)-mPerp*centerPoint(2);
    %sort data point base on x axis
    [x,idx] = sort(x);
    added = 4;
    
    
    if abs(mPerp) == inf
        yVec = centerPoint(1)-(ceil(vecLength/2+added)):centerPoint(1)+(ceil(vecLength/2+added));
        xVec = ones(size(yVec))*x(1);
    
    elseif abs(mPerp) <= 1
        xVec = centerPoint(2)-(ceil(vecLength+added)):centerPoint(2)+(ceil(vecLength+added));
        yVec = mPerp*xVec+b;   
        
    else
        yVec = centerPoint(1)-(ceil(vecLength+added)):centerPoint(1)+(ceil(vecLength+added));
        xVec = (yVec-b)/mPerp;

    end

    xVec = round(xVec);
    yVec = round(yVec);
       
    totalLength = (vecLength+added)*2;
    
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
