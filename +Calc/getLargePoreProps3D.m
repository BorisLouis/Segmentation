% The aim of this function is to calculate pore properties of e.g. polymer
% network images in 3D. 
% INPUT:
%       IMScaled: IM where the sampling in X,Y and Z is the same (in this case,
%                 pre-processed within mainPoreSizeCalc. This is crucial as
%                 otherwise the distance map will be odd as 1 px in XY
%                 would be different as 1 px in Z.
% OUTPUT:
%       pores3D: A structure containing the pore properties:
%               Inner and outer radius, connectivity, throats, position,
%               equivalent diameter,....
%      

function [pores3D] = getLargePoreProps3D(IM,pxSize)
    weights = [pxSize.XY,pxSize.XY,pxSize.Z];
    h = waitbar(0,'Searching for bubbles...');

    DMap = DistMap.calcWeightedDistMap(~IM,weights);
    pxDMap = bwdist(~IM);
    inRad3D = [];
    figure
    hold on
    view(3);
    
    while ~all(DMap(:)==0)
        [val] = max(DMap(:));
        
        idx = find(DMap==val);
        
        for i = 1:length(val)
            R = double(ceil(pxDMap(idx(i))));
            if length(inRad3D) >10
                %TODO test if value are still statistically significant
                
                if val < mean(inRad3D)-10*std(inRad3D)
                    disp('test')
                end
            end

            inRad3D = [inRad3D val];

            %Delete values in DMap in radius Val
            %create a sphere structural element 
            %R = ceil(val);
            SE = strel('sphere',R);
            
            %get position of the sphere centered around the maxima
            [y,x,z] = ind2sub(size(IM),idx(i));
            scatter(x,y,z,20,'filled')
            drawnow;
            
            xx = x-R:x+R;
            yy = y-R:y+R;
            zz = z-R:z+R;
            
            [xx,xid,yy,yid,zz,zid] = fixCoordinate(xx,yy,zz,size(IM));
            
            [X,Y,Z] = meshgrid(xx,yy,zz);
                       
            ind = sub2ind(size(IM),Y,X,Z);
            %create a volume with nothing
            tmpVol = zeros(size(IM));
            %replace the value by a sphere of 1 and 0s center around the
            %current maxima

            croppedSphere = SE.Neighborhood(yid,xid,zid);
            
            tmpVol(ind) = croppedSphere;
            
            DMap(logical(tmpVol)) = 0;
            
        end
        
        
        
 
        
        
    end
    pores3D.inRad = inRad3D;

end

function [X,X2,Y,Y2,Z,Z2] = fixCoordinate(X,Y,Z,dim)
    X2 = 1:length(X);
    Y2 = 1:length(Y);
    Z2 = 1:length(Z);
    
    
    %fix coordinates
    idx1 = X<1;
    X(idx1) = [];
    X2(idx1) = [];
    
    idx2 = Y<1;
    Y(idx2) = [];
    Y2(idx2) = [];
    
    idx3 = Z<1;
    Z(idx3) = [];
    Z2(idx3) = [];
    
    idx4 =(X>dim(2));
    X(idx4) = [];
    X2(idx4) = [];
   
    idx5 =(Y>dim(1));
    Y(idx5) = [];
    Y2(idx5) = [];
    
    idx6 =(Z>dim(3));
    Z(idx6) = [];
    Z2(idx6) = [];
    
   
 



end
