function [CDF, CCDF] = getCDF(dataV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUTHOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                            %
% Rafael Camacho Dejay (https://github.com/CamachoDejay)                  %                                                  %
% Website : Rafael Camacho Dejay: https://camachodejay.github.io/         %             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%
% In probability theory and statistics, the cumulative distribution
% function (CDF, also cumulative density function) of a real-valued random
% variable X, or just distribution function of X, evaluated at x, is the
% probability that X will take a value less than or equal to x.
dataV = double(dataV);
N = length(dataV);
d_sor = sort(dataV);
dt = diff(d_sor);
Ni = find(dt~=0);
Si = Ni./N;
Si = cat(1,0,Si,1);
xVal = unique(d_sor);
xVal = cat(1,xVal,max(xVal)+1);
% cleaning up edge effects
CDF.y = Si(2:end-1);
CDF.x = xVal(2:end-1);

% Sometimes, it is useful to study the opposite question and ask how often
% the random variable is above a particular level. This is called the
% complementary cumulative distribution function (CCDF) or simply the tail
% distribution or exceedance:
CCDF.y = 1-Si;
CCDF.x = xVal;
end

