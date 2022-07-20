%tmpPlot

figure

y = [Cold.porosity Hot.porosity;];
x = categorical({'Cold', 'Hot'});

bar(x,y);
ylim([0 1])
axis square
ylabel('Porosity');

%% Thickness1
minAll = min([ min(Cold.thicknessStats(:,1)) min(Hot.thicknessStats(:,1))]);
maxAll = max([ max(Cold.thicknessStats(:,1)) max(Hot.thicknessStats(:,1))]);

bins = minAll:300:maxAll;

[NCold] = histcounts(Cold.thicknessStats(:,1),bins);
[NHot] = histcounts(Hot.thicknessStats(:,1),bins);

figure
plot(bins(1:end-1),NCold/sum(NCold))
hold on
plot(bins(1:end-1),NHot/sum(NHot))

%% Thickness2
minAll = min([ min(Cold.thicknessStats(:,2)) min(Hot.thicknessStats(:,2))]);
maxAll = max([ max(Cold.thicknessStats(:,2)) max(Hot.thicknessStats(:,2))]);

bins = minAll:300:maxAll;

[NCold] = histcounts(Cold.thicknessStats(~isnan(Cold.thicknessStats(:,1)),2),bins);
[NHot] = histcounts(Hot.thicknessStats(:,2),bins);

figure
plot(bins(1:end-1),NCold/sum(NCold))
hold on
plot(bins(1:end-1),NHot/sum(NHot))

%% Straightness
minAll = min([ min(Cold.straightness) min(Hot.straightness)]);
maxAll = max([ max(Cold.straightness) max(Hot.straightness)]);

bins = minAll:0.05:maxAll;

[NCold] = histcounts(Cold.straightness(~isnan(Cold.straightness)),bins);
[NHot] = histcounts(Hot.straightness,bins);

figure
plot(bins(1:end-1),NCold/sum(NCold))
hold on
plot(bins(1:end-1),NHot/sum(NHot))
axis square
% 
% y = [0:0.1:1];
% x = repmat(median(Hot.straightness),1,length(y));
% hold on 
% plot(x,y,'r','Linewidth',1.5)
% y = [0:0.1:1];
% x = repmat(median(Cold.straightness),1,length(y));
% plot(x,y,'b','Linewidth',1.5)

%% Ratio
minAll = min([ min(Cold.ratio) min(Hot.ratio)]);
maxAll = max([ max(Cold.ratio) max(Hot.ratio)]);

bins = minAll:0.2:maxAll;

[NCold] = histcounts(Cold.ratio(~isnan(Cold.ratio)),bins);
[NHot] = histcounts(Hot.ratio,bins);

figure
plot(bins(1:end-1),NCold/sum(NCold))
hold on
plot(bins(1:end-1),NHot/sum(NHot))
axis square

% y = [0:0.1:1];
% x = repmat(median(Hot.straightness),1,length(y));
% hold on 
% plot(x,y,'r','Linewidth',1.5)
% y = [0:0.1:1];
% x = repmat(median(Cold.straightness),1,length(y));
% plot(x,y,'b','Linewidth',1.5)





