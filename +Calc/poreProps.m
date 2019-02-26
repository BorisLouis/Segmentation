% The aim of this function is to calculate the pore size and other properties 
% of pores based on a pre-segmented image.

%INPUT:
%       IM: A binary image or binary volume
%       IMScaled: IM where the sampling in X,Y and Z is the same (processed
%                 within mainPoreSizeCalc(only needed for 3D)
%       dim: '2D', '3D' and 'both' user option to run the analysis only in
%       a specific dimension (2D analyze z-slice independently, 3D analyze
%       considering the data as a coherent volume. 'both' will run the two
%       types
%             
%OUTPUT:
%       pores2D: A strutcture containing different poreProperties:
%                inner and external radii, throats, connectivity, equiv
%                diameter, center position
%                return empty if 2D analysis is not run
%       pores3D: A structure containing the 3D equivalent of the pores2D
%                properties
%                return empty if 2D analysis is not run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boris Louis (https://github.com/BorisLouis)                             %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %
% Susana Rocha (https://github.com/SusanaRocha)                           %                                                                        %
% Website : Boris Louis: https://borislouis.github.io/                    %
%           Rafael Camacho Dejay: https://camachodejay.github.io/         %                               
%           Susana Rocha: https://susanarocha.github.io/                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [pores2D,pores3D] = poreProps(IM,IMScaled,dim)
%check number of input provided
switch nargin
    case 1
        dim = '2D'; %if IM is not scaled 3D should not be used
    case 2
        dim = 'both';
    case 3
    otherwise
        error('Wrong number of input argument, expect [1 3]');
end
% act depending on the dimension requested by the user
switch dim
    case '2D'
        an2D = true;
        an3D = false;
        pores3D =[];
    case '3D'
        an2D = false;
        an3D = true;
        pores2D = [];
    case 'both'
        an2D = true;
        an3D = true;
    otherwise
        error('Please provided the type of analysis you want to perform (2D, 3D or both)');
end

nFrames = size(IM,3);

if an2D
    % create variables for storing
    all_areas = [];
    all_inRad = [];
    all_extRad = [];
    all_throats = [];
    all_connect = [];
    all_ctr_ext = []; 
    diameter = [];
    
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
        ws = bwareaopen(ws,2);%remove pores smaller than 2px
        ws = bwlabel(ws);
%      %   figure, imagesc(ws);

        %% calculate data
        %Calculate pore properties in 2D
        stats = regionprops(logical(ws),'PixelIdxList','Area', 'Image', 'Centroid','EquivDiameter', 'BoundingBox'); %ADDED BOUNDINGBOX SUSANA
        %memory preallocation
        extRad = zeros(1, max(ws(:)));
        inRad = zeros(1, max(ws(:)));
        area = zeros(1, max(ws(:)));
        ctr = zeros(max(ws(:)),2); %center, for connectivity
        ctr_ext = zeros(max(ws(:)),3); %center of circle/sphere and frame
        NRconn = zeros(1, max(ws(:))); %nr of connected pores
        %process pore individually to extract additional parameters
        for i = 1:max(ws(:))
            extRad(i) = max(D(ws==i));
            inRad(i) = max(max(bwdist(~stats(i).Image)));
            area(i) = stats(i).Area;
            ctr(i, :) = stats(i).Centroid;
            %Calculate connectivity by expanding the coordinate of a pores
            %and checking overlap with others
            [I,J] = ind2sub(size(im),stats(i).PixelIdxList);
            coord = [I-ctr(i,2) J-ctr(i,1)];
            scaleFactor = (abs(coord)+2)./abs(coord);
            %Calculation based on coordinate because much faster to process
            %than using imdilate on the full image/volume
            dilCoord = coord;
            dilCoord(dilCoord(:,1)~=0,1) = coord(coord(:,1)~=0,1).*scaleFactor(coord(:,1)~=0,1);
            dilCoord(dilCoord(:,2)~=0,2) = coord(coord(:,2)~=0,2).*scaleFactor(coord(:,2)~=0,2);
    
            dilCoord  =[dilCoord(:,1)+ctr(i,2) dilCoord(:,2)+ctr(i,1)];
            dilCoord(dilCoord<1) = 1;
            dilCoord(dilCoord>size(im,1)) = size(im,1);
            idx = sub2ind(size(im),dilCoord(:,1),dilCoord(:,2));
            idx = round(idx);
            
            conn_pores = ws(idx);
            conn_pores(conn_pores==0) = [];

            IDconn = unique(conn_pores);
            NRconn(i) = length(IDconn);
            
            % calculate position of max distance (use extRad) for plotting
            % pores on 3D model
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
end

if an3D
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
     
    %individual treatment of the pores
    for i = 1:size(stats,1)
        extRad3D(i) = max(D(stats(i,:).VoxelIdxList{1}));
        inRad3D(i) = max(max(max(bwdist(~stats(i,:).Image{1}))));
        vol(i) = stats(i,:).Volume;
        ctr3D(i, :) = stats(i,:).Centroid;
        
        [I,J,K] = ind2sub(size(IM),stats(i,:).VoxelIdxList{1});
        coord = [I-ctr3D(i,2) J-ctr3D(i,1) K-ctr3D(i,3)];
        scaleFactor = (abs(coord)+2)./abs(coord);
        dilCoord = coord;
        dilCoord(dilCoord(:,1)~=0,1) = coord(coord(:,1)~=0,1).*scaleFactor(coord(:,1)~=0,1);
        dilCoord(dilCoord(:,2)~=0,2) = coord(coord(:,2)~=0,2).*scaleFactor(coord(:,2)~=0,2);
        dilCoord(dilCoord(:,3)~=0,3) = coord(coord(:,3)~=0,3).*scaleFactor(coord(:,3)~=0,3);
        
        dilCoord  =[dilCoord(:,1)+ctr3D(i,2) dilCoord(:,2)+ctr3D(i,1) dilCoord(:,3)+ctr3D(i,3)];
        dilCoord(dilCoord<1) = 1;
        dilCoord(dilCoord(:,1)>size(IM,1),1) = size(IM,1);
        dilCoord(dilCoord(:,2)>size(IM,2),2) = size(IM,2);
        dilCoord(dilCoord(:,3)>size(IM,3),3) = size(IM,3);
        
        idx = sub2ind(size(IM),dilCoord(:,1),dilCoord(:,2),dilCoord(:,3));
        idx = round(idx);
        
        conn_pores = ws(idx);
        conn_pores(conn_pores==0) = [];

        IDconn = unique(conn_pores);
        NRconn(i) = length(IDconn);
        
        % calculate position of max distance (use extRad) ADDED SUSANA
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
        
%         if any(ctr_ext3D(i,:)>size(IM,1))
%             disp('HAAAAA');
%         end

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
    pores3D.ctr_ext = ctr_ext3D; %ADDED SUSANA
    h = waitbar(1,h,'Done');

    pause(1);

    close(h);
end


