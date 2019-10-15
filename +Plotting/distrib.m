function distrib(data,centers,log)
    %get bin
    d = diff(centers)/2;
    edges = [centers(1)-d(1), centers(1:end-1)+d, centers(end)+d(end)];
    %create histogram data
    [N, ~]=histcounts(data, edges);
    Nf = {cat(1,centers,N)'};
    %plot distribution
    Plotting.distributionPlot(Nf,'histOpt',0, 'showMM',0);
    pp=findobj(gca, 'type', 'patch');
    set(pp(1), 'FaceColor', [0 0 1]);    
    hold on
    
    h = boxplot(data, 'OutlierSize',0.1, 'Colors', [0.2 0.2 0.2; 0.2 0.2 0.2], 'Widths',0.9);
    
    set(h([1:5,7],:),'Visible','off')
    set(h(6, :), 'LineWidth', 2)
    if log
        set(gca, 'YScale', 'log')
    end
 
   

end