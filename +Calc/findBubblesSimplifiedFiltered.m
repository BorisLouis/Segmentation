function [bubble_radii, bubble_coord] = findBubblesSimplifiedFiltered(img)
% find_bubbles determines the location and radii of 2D or 3D bubbles in
% the binary image img.
%
% 19|03|2013 Stefan M�nster
EDM = bwdist(img); % This calculation of the EDM of the fluid phase assumes
%that fluid pixels are 0 and solid pixels are 1;
smoothed_EDM = imfilter(EDM,fspecial('gaussian',5,1)); % this smoothing
%suppresses bubbles of similar size in close proximity and can be omitted
local_maxima = imregionalmax(smoothed_EDM); % determines the local maxima,
%which are the center locations of the bubbles

if(length(size(img))<3) % translates local maxima into x,y,z coordinates
    
 [bubble_coord(:,1),bubble_coord(:,2)]=...
     ind2sub(size(img),find(local_maxima));
 
else
    
 [bubble_coord(:,1), bubble_coord(:,2),bubble_coord(:,3)]=...
     ind2sub(size(img),find(local_maxima));
 
end

%Filtering the data

bubble_radii=EDM(local_maxima); % determines the radii of the bubbles from
%the EDM values at the local maxima
end