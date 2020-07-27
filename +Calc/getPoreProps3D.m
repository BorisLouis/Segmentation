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

function [pores3D] = getPoreProps3D(IM,pxSize)
%% Performing analysis in 3D
%The calculation here are the same than the one in 2D but some
%adaptation was required to perform it in 3D
%Gaussian filtering
weights = [pxSize.XY,pxSize.XY,pxSize.Z];
h = waitbar(0,'3D Gaussian filtering...');

DMap = DistMap.calcWeightedDistMap(~IM,weights);
r = 3;
Ds = imgaussfilt3(DMap, [r r r*pxSize.XY/pxSize.Z]);

%previous code
% D = bwdist(~IMScaled);

% Ds = imgaussfilt3(D, [r r r]);

%Local maxima and merging
waitbar(0.2,h,'Getting Local maxima...');
LocMax = imregionalmax(Ds); %find local maxima
r=10;
se = strel('sphere',r);
LocMax2 = imdilate(LocMax,se);
r=8;
se = strel('sphere',r);
LocMax3 = imerode(LocMax2,se);
LocMax3(~IM) = 0;

%release some memory
clear LocMax LocMax2

waitbar(0.4,h,'Performing watershed...');
imD = -Ds;
%clear memory
clear Ds

imD(~IM) = Inf;
imD(LocMax3) = min(imD(:));
ws = watershed(imD);
ws(~IM) = 0;
ws = bwareaopen(logical(ws),4);%remove small pores
ws = bwlabeln(logical(ws));

%clear memory
clear imD LocMax3

h = waitbar(0.6,h,'Getting pore 3D properties...');
%Calculation of pore properties in 3D:
stats = regionprops3(ws, 'VoxelIdxList', 'Volume', 'Image', 'Centroid', 'BoundingBox','EquivDiameter'); %ADDED SUSANA
%Memory preallocation
extRad3D = zeros(size(stats,1),1);
inRad3D  = zeros(size(stats,1),1);
vol      = zeros(size(stats,1),1);
NRconn   = zeros(size(stats,1),1);
ctr3D    = zeros(size(stats,1),3);
ctr_ext3D  = zeros(size(stats,1),3);
IDconn = cell(1,max(ws(:)));

%individual treatment of the pores
for i = 1:size(stats,1)
    extRad3D(i) = max(DMap(stats(i,:).VoxelIdxList{1}));
    edmIm = DistMap.calcWeightedDistMap(~stats(i,:).Image{1},weights);
    inRad3D(i) = max(edmIm(:));
   % inRad3D(i) = max(max(max(bwdist(~stats(i,:).Image{1}))));
    vol(i) = stats(i,:).Volume;
    ctr3D(i, :) = stats(i,:).Centroid;
    %Calculate connectivity by expanding the coordinate of a pores
    %and checking overlap with others
    
    poreImage = stats(i,:).Image{1};
    dilated   = zeros(size(poreImage,1)+6,size(poreImage,2)+6,...
        size(poreImage,3)+6);
    
    dilated(4:end-3,4:end-3,4:end-3) = poreImage;
    r=2;
    se = strel('sphere',r);
    dilated = imdilate(dilated,se);
    Coord = find(dilated==1);
    
    [I,J,K] = ind2sub(size(dilated),Coord);
    dilCoord = [I(:),J(:),K(:)];
    ctrTmp = mean(dilCoord,1);
    
    dilCoord = dilCoord-ctrTmp + ctr3D(i,[2 1 3]);

    %clean borders
    dilCoord(dilCoord<1) = 1;
    dilCoord(dilCoord(:,1)>size(DMap,1),1) = size(DMap,1);
    dilCoord(dilCoord(:,2)>size(DMap,2),2) = size(DMap,2);
    dilCoord(dilCoord(:,3)>size(DMap,3),3) = size(DMap,3);
    %round for indexing
    dilCoord = round(dilCoord);
    
    idx = sub2ind(size(DMap),dilCoord(:,1),dilCoord(:,2),dilCoord(:,3));
    idx = round(idx);
    %get the connected pores
    conn_pores = ws(idx);
    conn_pores(conn_pores==0) = [];
    %find which pore are connected and their numbers
    connID    = unique(conn_pores);
    IDconn{i} = connID;
    NRconn(i) = length(connID);
    
    % calculate position of max distance (use extRad)
    x = round(stats.BoundingBox(i,1));
    y = round(stats.BoundingBox(i,2));
    z = round(stats.BoundingBox(i,3));
    wx = stats.BoundingBox(i,4);
    wy = stats.BoundingBox(i,5);
    wz = stats.BoundingBox(i,6);
    
    sel_im = DMap(y:y+wy-1,x:x+wx-1, z:z+wz-1);
    clear ind
    idx=find(sel_im==max(sel_im(:)));
    [ind(:,2),ind(:,1),ind(:,3)] = ind2sub(size(sel_im),idx);
    ctr_ext3D(i,:)=round(stats.BoundingBox(i,1:3))+mean(ind,1)-1;
    
end
%clear memory
clear stats x y z wx wy wz
%need to fi all of it below
inRad3D(inRad3D == Inf) = NaN;

ws = logical(ws); % used to calculate overlap reagions
mask_throats = IM-ws; % mask for the throats
throats = table2array(regionprops3(logical(mask_throats),DMap,'MaxIntensity'));

h = waitbar(0.9,h,'Storing data');
pores3D.vol = vol;
pores3D.extRad = extRad3D;
pores3D.inRad = inRad3D;
pores3D.throats = throats;
pores3D.connect = NRconn;
pores3D.ctr_ext = ctr_ext3D;
pores3D.connID  = IDconn;

h = waitbar(1,h,'Done');

pause(1);

close(h);
end