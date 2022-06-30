%% Simplified bubble analysis
% This code was taken from the supplementary information of the following
% paper:
% Munster & Fabry : A Simplified implementation of the bubble analysis of
% biopolymer network
% the paper can be read here : 
% https://www.cell.com/biophysj/fulltext/S0006-3495(13)00572-9

function [bubble_radii, bubble_coord] = findBubblesSimplified(img,pxSize)
% find_bubbles determines the location and radii of 2D or 3D bubbles in
% the binary image img.
%
% 19|03|2013 Stefan Münster
%EDM = bwdist(img);% This calculation of the EDM of the fluid phase assumes
%that fluid pixels are 0 and solid pixels are 1;
weight(1) = pxSize.XY;
weight(2) = pxSize.XY;
weight(3) = pxSize.Z;
EDM = DistMap.calcWeightedDistMap(img,weight);

smoothed_EDM = imfilter(EDM,fspecial('gaussian',10,1)); % this smoothing
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
bubble_radii=EDM(local_maxima); % determines the radii of the bubbles from
%the EDM values at the local maxima
end