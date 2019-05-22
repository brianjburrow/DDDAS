function niceErrorbar(x, y, std)
    dy = 3*std;
    h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 .9 .9],'linestyle','none');
    set(h,'facealpha',.85)
    %line(x,y, 'LineWidth', 2.0)
end