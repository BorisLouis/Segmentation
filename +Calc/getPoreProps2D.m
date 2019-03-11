function [pores2D] = getPoreProps2D(IM)
% create variables for storing
    all_areas = [];
    all_inRad = [];
    all_extRad = [];
    all_throats = [];
    all_connect = [];
    all_ctr_ext = [];
    all_connectID = [];
    diameter = [];
    nFrames = size(IM,3);
    h=waitbar(0, 'Calculating watershed / analysing pore diameter in 2D...');
    for fr = 1 : nFrames
        im = IM(:,:,fr); %im is the binary image from analysis

      % calculate distance matrix
        D = bwdist(~im);

      % smooth distance matrix with gaussian filter
        r = 3;
        Ds = imgaussfilt(D, r);

        % merge local maxima
        LocMax = imregionalmax(Ds); %find local maxima
        LocMax2 = imdilate(LocMax, strel('disk',10));%merge local maxima
        LocMax3 = bwmorph(LocMax2,'thin',8);
        LocMax3(~im) = 0;
    %    figure, imshowpair(LocMax3,LocMax2,'falsecolor')

        %% perform watershed
        imD = -Ds;
        imD(~im) = Inf;
        imD(LocMax3) = min(imD(:));
        ws = watershed(imD);
        ws(~im) = 0;
        ws = bwareaopen(logical(ws),2);%remove pores smaller than 2px
        ws = bwlabel(logical(ws));
%      %   figure, imagesc(ws);

        %% calculate data
        %Calculate pore properties in 2D
        stats = regionprops(ws,'PixelIdxList','Area', 'Image', 'Centroid','EquivDiameter', 'BoundingBox'); %ADDED BOUNDINGBOX SUSANA
        %memory preallocation
        extRad = zeros(1, max(ws(:)));
        inRad = zeros(1, max(ws(:)));
        area = zeros(1, max(ws(:)));
        ctr = zeros(max(ws(:)),2); %center, for connectivity
        ctr_ext = zeros(max(ws(:)),3); %center of circle/sphere and frame
        IDconn = cell(1,max(ws(:)));
        NRconn = zeros(1, max(ws(:))); %nr of connected pores
        %process pore individually to extract additional parameters
        for i = 1:max(ws(:))
            extRad(i) = max(D(ws==i));
            inRad(i) = max(max(bwdist(~stats(i).Image)));
            area(i) = stats(i).Area;
            ctr(i, :) = stats(i).Centroid;
            
            
            %Calculate connectivity by expanding the coordinate of a pores
            %and checking overlap with others
            poreImage = stats(i,:).Image;
            dilated   = zeros(size(poreImage,1)+6,size(poreImage,2)+6);
    
            dilated(4:end-3,4:end-3) = poreImage;
            r=2;
            se = strel('disk',r);
            dilated = imdilate(dilated,se);
            Coord = find(dilated==1);

            [I,J] = ind2sub(size(dilated),Coord);
            dilCoord = [I(:),J(:)];
            ctrTmp = mean(dilCoord,1);

            dilCoord = dilCoord-ctrTmp + ctr(i,[2 1]);

            %clean borders
            dilCoord(dilCoord<1) = 1;
            dilCoord(dilCoord(:,1)>size(im,1),1) = size(im,1);
            dilCoord(dilCoord(:,2)>size(im,2),2) = size(im,2);
            
            %round for indexing
            dilCoord = round(dilCoord);

            idx = sub2ind(size(im),dilCoord(:,1),dilCoord(:,2));
            idx = round(idx);
            %%%%%%%
%             [I,J] = ind2sub(size(im),stats(i).PixelIdxList);
%             coord = [I-ctr(i,2) J-ctr(i,1)];
%             scaleFactor = (abs(coord)+2)./abs(coord);
%             %Calculation based on coordinate because much faster to process
%             %than using imdilate on the full image/volume
%             dilCoord = coord;
%             dilCoord(dilCoord(:,1)~=0,1) = coord(coord(:,1)~=0,1).*scaleFactor(coord(:,1)~=0,1);
%             dilCoord(dilCoord(:,2)~=0,2) = coord(coord(:,2)~=0,2).*scaleFactor(coord(:,2)~=0,2);
%     
%             dilCoord  =[dilCoord(:,1)+ctr(i,2) dilCoord(:,2)+ctr(i,1)];
%             dilCoord(dilCoord<1) = 1;
%             dilCoord(dilCoord>size(im,1)) = size(im,1);
%             idx = sub2ind(size(im),dilCoord(:,1),dilCoord(:,2));
%             idx = round(idx);
            
            conn_pores = ws(idx);
            conn_pores(conn_pores==0) = [];
            connID = unique(conn_pores);
            frame = ones(size(connID))*fr;
            IDconn{i} = [connID frame];
            NRconn(i) = length(connID);
            
            % calculate position of max distance (use extRad) for plotting
            % pores 
            x = round(stats(i).BoundingBox(1));
            y = round(stats(i).BoundingBox(2));
            wx = stats(i).BoundingBox(3);
            wy = stats(i).BoundingBox(4);
            sel_im = D(y:y+wy-1,x:x+wx-1);
            clear ind
            [ind(:,2), ind(:,1)]=find(sel_im==max(sel_im(:)));
            ctr_ext(i,:)=[stats(i).BoundingBox(1:2)+mean(ind,1) fr];

        end
        %clean data
        inRad(inRad==Inf) = NaN;
        inRad(isnan(inRad)) = [];
        %calculation of the throats
        ws1 = logical(ws); % used to calculate overlap reagions
        mask_throats = im-ws1; % mask for the throats
        throats = struct2array(regionprops(logical(mask_throats),D,'MaxIntensity'));
        
        %save data for each frame
        all_areas   = [all_areas area];
        all_inRad   = [all_inRad inRad];
        all_extRad  = [all_extRad extRad];
        all_throats =  [all_throats throats];
        all_connect = [all_connect NRconn]; 
        all_connectID = {all_connectID IDconn};
        all_ctr_ext = [all_ctr_ext; ctr_ext];  % added SUSANA
        diameter = [diameter cell2mat({stats.EquivDiameter})];
        %solidity = [solidity cell2mat({stats.Solidity})];
        
        waitbar(fr/nFrames, h);
    end
    close (h)
    %storing as function output
    pores2D.area = all_areas;
    pores2D.inRad = all_inRad;
    pores2D.extRad = all_extRad;
    pores2D.throats = all_throats;
    pores2D.connect = all_connect;
    pores2D.diameter = diameter;
    pores2D.ctr_ext = all_ctr_ext;
    pores2D.connID = all_connectID;
end