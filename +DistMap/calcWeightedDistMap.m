function [EDM] = calcWeightedDistMap(mask,weight)
  
    mask = logical(mask);
    
    [test,idx] = bwdist(mask);
   
    % following step are in a loop to be more memory efficient:
        EDM = zeros(size(mask));
        for i = 1 : size(mask,3)

            [row,col,z] = ind2sub(size(mask),idx(:,:,i));
            [X,Y] = meshgrid(1:size(row,2),1:size(col,1));

            EDM(:,:,i) = sqrt(((X-col)*weight(1)).^2 + ...
                   ((Y-row)*weight(2)).^2 + ...
                   ((i-z)*weight(3)).^2);
        end
end