function plot3dhist(samples, numSlice, mu, figNum, standNorm)
    for idx = 1:numSlice
        disp(idx)
        hfig = figure(12524);
        hist(samples(:,idx), 25, 'black');
        h = get(gca,'Children'); 
        x = h.Vertices(:,1);
        z = h.Vertices(:,2)./max(h.Vertices(:,2));
        y = mu(idx)*ones(size(x));
        figure(figNum)
        patch(x,y,z, [0.6350, 0.0780, 0.1840]); hold on; view([56.1, 82.8])
    end
    figure(figNum)
    zlabel("Posterior Hist.")
    ylabel("Measurement")
    xlabel("Unknown State")
    grid on
    %AxesH = axes('Ylim', [0, 4], 'YTick', 0:0.1:4, 'NextPlot', 'add');
    %view([-55.9000, 68.4])
     tm = [min(min(samples)), max(max(samples)), max(max(samples)), min(min(samples))];
     ty = [min(mu), min(mu), max(mu), max(mu)];
     tz = [0,0,0,0];
     p = patch(tm, ty, tz, 'black');
     set(p,'facealpha',0.1)
     set(p,'edgealpha',0.1)
     if standNorm
         xlim([-4, 4])
     else
        xlim([min(min(samples)), max(max(samples))])
     end
     
end