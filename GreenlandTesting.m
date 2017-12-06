setworkspace()
G = glmalpha('greenland',90);
lmcosi = grabanalpha(G,N,L,1,2)
%%
% k = 1(1,0) 2(1,1) 5(1,-1) 7(2,1)
for k = 26:30
    figure(k)
    N = 8;L=90;
    % plot a Greenland Eigentaper
    grabanalpha(G,N,L,k,2)
    hold on
    
    % improvise a grid
    lon_degree = 6;
    lat_degree = lon_degree/2;
    sphere(360/lon_degree);
    
    % focus on Greenland
    view(45,90)
    map_r = 0.4;
    axis([-0 map_r -map_r 0])
    % [-0 map_r -map_r 0]
    % [-1 1 -1 1]
    
    % turn on axis, turn off ticks and title
    axis on
    title('')
    set(gca,'xtick',[],'ytick',[])
    
    % longitude ticks
    degy = (-7:0)*lon_degree;
    degx = (0:7)*lon_degree;
    y = [map_r*tan(deg2rad(degy))];
    x = [map_r*tan(deg2rad(degx))];
    xlab = strcat('$',num2str(-90+degx'),'^{\circ}$');
    ylab = strcat('$',num2str(degy'),'^{\circ}$');
    for i = 1:length(degy)
        %ylabels
        text(map_r+0.02,y(i),0,ylab(i,:),'interpreter','latex','horizontalalignment','center')
        text(x(i),-(map_r+0.02),0,xlab(i,:),'interpreter','latex','horizontalalignment','center')
    end
    
    % latitude ticks
    degz = (0:7).*lat_degree;
    x = sin(deg2rad(degz));
    z = cos(deg2rad(degz));
    for i = 1:length(degz)
        text(x(i),0+0.01,0,strcat(num2str(degz(i)),'$^{\circ}$'),'interpreter','latex')
    end
    
    % axis labeling
    text(map_r+0.055,-0.5*map_r,1,'Longitude','horizontalalignment','center','rotation',45);
    text(0.5*map_r,0.05,1,'Colatitude','horizontalalignment','center','rotation',-45);
    hold off
end