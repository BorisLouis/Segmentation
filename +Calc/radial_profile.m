function [Tics,Average]=radial_profile(data,radial_step)
%main axii cpecified:
x=(1:size(data,2))-(size(data,2)+1)/2;
y=(1:size(data,1))-(size(data,1)+1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);
end