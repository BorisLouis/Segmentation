function [pores3D] = getPoreProps3D(IMScaled) 

%% Performing analysis in 3D
    %The calculation here are the same than the one in 2D but some
    %adaptation was required to perform it in 3D
    %Gaussian filtering
    h = waitbar(0,'3D Gaussian filtering...');
    D = bwdist(~IMScaled);
    r = 3;
    Ds = imgaussfilt3(D, [r r r]);
    
    %Local maxima and merging
    waitbar(0.2,h,'Getting Local maxima...');
    LocMax = imregionalmax(Ds); %find local maxima
    r=10;
    se = strel('sphere',r);
    LocMax2 = imdilate(LocMax,se);
    r=8;
    se = strel('sphere',r);
    LocMax3 = imerode(LocMax2,se); 
    LocMax3(~IMScaled) = 0;
    waitbar(0.4,h,'Performing watershed...');
    imD = -Ds;
    imD(~IMScaled) = Inf;
    imD(LocMax3) = min(imD(:));
    ws = watershed(imD);
    ws(~IMScaled) = 0;
    ws = bwareaopen(ws,4);%remove small pores
    ws = bwlabeln(ws);
    %figure, imagesc(ws(:,:,6));
    h = waitbar(0.6,h,'Getting pore 3D properties...');
    %Calculation of pore properties in 3D:
    stats = regionprops3(logical(ws), 'voxelIdxList', 'Volume', 'Image', 'Centroid','EquivDiameter', 'BoundingBox'); %ADDED SUSANA
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
        extRad3D(i) = max(D(stats(i,:).VoxelIdxList{1}));
        inRad3D(i) = max(max(max(bwdist(~stats(i,:).Image{1}))));
        vol(i) = stats(i,:).Volume;
        ctr3D(i, :) = stats(i,:).Centroid;
         %Calculate connectivity by expanding the coordinate of a pores
            %and checking overlap with others
        [I,J,K] = ind2sub(size(IMScaled),stats(i,:).VoxelIdxList{1});
        coord = [I-ctr3D(i,2) J-ctr3D(i,1) K-ctr3D(i,3)];
        %Scale factor is define for each pore to get it about 2 pixel
        %bigger than it really is (hence connecting it to other pores as if
        %we had not performed watershed)
        scaleFactor = (abs(coord)+2)./abs(coord);
        dilCoord = coord;
        dilCoord(dilCoord(:,1)~=0,1) = coord(coord(:,1)~=0,1).*scaleFactor(coord(:,1)~=0,1);
        dilCoord(dilCoord(:,2)~=0,2) = coord(coord(:,2)~=0,2).*scaleFactor(coord(:,2)~=0,2);
        dilCoord(dilCoord(:,3)~=0,3) = coord(coord(:,3)~=0,3).*scaleFactor(coord(:,3)~=0,3);
        %Calculation based on coordinate because much faster to process
        %than using imdilate on the full image/volume 
        dilCoord  =[dilCoord(:,1)+ctr3D(i,2) dilCoord(:,2)+ctr3D(i,1) dilCoord(:,3)+ctr3D(i,3)];
        dilCoord(dilCoord<1) = 1;
        dilCoord(dilCoord(:,1)>size(IMScaled,1),1) = size(IMScaled,1);
        dilCoord(dilCoord(:,2)>size(IMScaled,2),2) = size(IMScaled,2);
        dilCoord(dilCoord(:,3)>size(IMScaled,3),3) = size(IMScaled,3);
        
        idx = sub2ind(size(IMScaled),dilCoord(:,1),dilCoord(:,2),dilCoord(:,3));
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
        sel_im = D(y:y+wy-1,x:x+wx-1, z:z+wz-1);
        clear ind
        idx=find(sel_im==max(sel_im(:)));
        [ind(:,2),ind(:,1),ind(:,3)] = ind2sub(size(sel_im),idx);
        ctr_ext3D(i,:)=round(stats.BoundingBox(i,1:3))+mean(ind,1)-1;
    end

    %need to fi all of it below
    inRad3D(inRad3D==Inf) = NaN;
    
    ws1 = logical(ws); % used to calculate overlap reagions
    mask_throats = IMScaled-ws1; % mask for the throats
    throats = table2array(regionprops3(logical(mask_throats),D,'MaxIntensity'));
    h = waitbar(0.9,h,'Storing data');
    pores3D.vol = vol;
    pores3D.extRad = extRad3D;
    pores3D.inRad = inRad3D;
    pores3D.throats = throats;
    pores3D.connect = NRconn;
    pores3D.diameter = stats.EquivDiameter;
    pores3D.ctr_ext = ctr_ext3D;
    pores3D.connID  = IDconn;
    h = waitbar(1,h,'Done');

    pause(1);

    close(h);
end