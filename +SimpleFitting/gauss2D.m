function [G] = gauss2D(pos, sig, xid,yid,maxCount)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here 
switch nargin
    case 4
        A = 100;
    case 5
        A = maxCount;
    otherwise
        error('Not enough input');
end

[x,y] = meshgrid(xid,yid);
sigX = sig(1);
sigY = sig(2);
x0 = pos(1);
y0 = pos(2);

xPart = ((x-x0).^2) ./ (2*sigX^2);
yPart = ((y-y0).^2) ./ (2*sigY^2);


G = A.*exp( -(xPart + yPart));

end