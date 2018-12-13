function latlonaxis( xlim, ylim, thespacelim, ntix )
%LATLONAXIS command for lat/lon axis lables

xtix = linspace(xlim(1),xlim(2),ntix);
xtiklab = linspace(thespacelim(1,1),thespacelim(1,2),ntix);
ytix = linspace(ylim(1),ylim(2),ntix);
ytiklab = linspace(thespacelim(2,2),thespacelim(2,1),ntix);

xtiklab = strcat(string(xtiklab),'\circ');
ytiklab = strcat(string(ytiklab),'\circ');

latlonax =  strcat("set(gca,",...
    "'xtick',xtix,'xticklabel',xtiklab,",...
    "'ytick',ytix,'yticklabel',ytiklab,",...
    "'xgrid','on','ygrid','on',",...
    "'xminorgrid','on','yminorgrid','on',",...
    "'minorgridcolor','w','gridcolor','w')");
eval(latlonax);
xlabel('longitude');
ylabel('latitude');
end

