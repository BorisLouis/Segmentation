function [ mov ] = getframes( path2file, frames )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %                                                  %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%   Detailed explanation goes here

assert(min(size(frames))==1,'frames must be a vector of positive integers')

f_n = length(frames);

tObj = Tiff(path2file,'r');
l    = tObj.getTag(256);
w    = tObj.getTag(257);
tObj.setDirectory(frames(1));

im1  = tObj.read;
nClass = class(im1);
mov = zeros(w,l,f_n,nClass);
convert = false;

if length(size(im1))~=2
    warning('Colored images received, conversion to grayscale is performed')
    convert = true;
end

for i = 1:f_n
    f_i = frames(i);
    tObj.setDirectory(f_i)
    movTmp = tObj.read;  
    if convert
        movTmp = rgb2gray(movTmp);
    end
    mov(:,:,i) = movTmp;    
end
tObj.close


end

