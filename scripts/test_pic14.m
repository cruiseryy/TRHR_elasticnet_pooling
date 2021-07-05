function test_pic14(map)
% figure('Color','white');
% set(gcf,'outerposition',get(0,'screensize'));
cmax = max(max(map));
cmin = min(min(map));
chalf = max((cmax - cmin)/2,1e-8);
ccenter = 0.5*(cmax+cmin);
set(gca,'CLim',[ccenter-chalf, ccenter+chalf])
% set(gca, 'CLim', [-2,1.5]);  
% map = zeros(360,80);
ax=axesm('MapProjection','eqdcylin','MapLatLimit',[-90 90],'MapLonLimit',[0 360]);
setm(ax,'MlabelLocation', 30,'PlabelLocation',[-90 -60 -30 0 30 60 90],'MeridianLabel','on','ParallelLabel','on','MLabelParallel','north');
% setm(gca,'grid','on')
R = georasterref('RasterSize', size(map'),'Latlim', [-90 90], 'Lonlim', [0 360]);
geoshow(map',R, 'DisplayType','texturemap')
% geoshow(lat1,lon1,'DisplayType','line')
load('coast');
geoshow(lat,long,'DisplayType','line','Color','black','LineWidth',1);

% mincolor = [0 0.45 0.74];
% maxcolor = [0.6 0.2 0];
% middlecolor = [1 1 1];
% customcmap1 = [linspace(mincolor(1),middlecolor(1),6)',...
%             linspace(mincolor(2),middlecolor(2),6)',...
%             linspace(mincolor(3),middlecolor(3),6)'];
% customcmap2 = [linspace(middlecolor(1),maxcolor(1),6)',...
%             linspace(middlecolor(2),maxcolor(2),6)',...
%             linspace(middlecolor(3),maxcolor(3),6)'];
% customcmap = [customcmap1(1:5,:); customcmap2];
% colormap(customcmap)
% colorbar

% mincolor = [0 0.75 0.75];
% maxcolor = [0.87 0.49 0];
% middlecolor = [1 1 1];
% customcmap1 = [linspace(mincolor(1),middlecolor(1),10)',...
%             linspace(mincolor(2),middlecolor(2),10)',...
%             linspace(mincolor(3),middlecolor(3),10)'];
% customcmap2 = [linspace(middlecolor(1),maxcolor(1),10)',...
%             linspace(middlecolor(2),maxcolor(2),10)',...
%             linspace(middlecolor(3),maxcolor(3),10)'];
% customcmap = [customcmap1(1:9,:); customcmap2];
% colormap(customcmap)
mincolor = [0 0.75 0.75];
maxcolor = [0.87 0.49 0];
middlecolor = [1 1 1];
customcmap1 = [linspace(mincolor(1),middlecolor(1),10)',...
            linspace(mincolor(2),middlecolor(2),10)',...
            linspace(mincolor(3),middlecolor(3),10)'];
customcmap2 = [linspace(middlecolor(1),maxcolor(1),10)',...
            linspace(middlecolor(2),maxcolor(2),10)',...
            linspace(middlecolor(3),maxcolor(3),10)'];
customcmap = [customcmap1(1:9,:); customcmap2];
colormap(customcmap)
colorbar




% geoshow(data(:,2),data(:,1),'DisplayType','line','Color','black','LineWidth',1);
geoshow([-18 -18 -42 -42 -18],[180 200 200 180 180],'DisplayType','Line','color','k','linewidth',1.5);
% geoshow([-50 -50],[180 360 ],'DisplayType','Line','color','k','linewidth',2);
% geoshow([-50 -50],[0 180],'DisplayType','Line','color','k','linewidth',2);
% geoshow([50 50],[180 360 ],'DisplayType','Line','color','k','linewidth',2);
% geoshow([50 50],[0 180],'DisplayType','Line','color','k','linewidth',2);
% geoshow([-50 50],[359.99 359.99],'DisplayType','Line','color','k','linewidth',1);
% geoshow([-5 -5],[1 359.99],'DisplayType','Line','color','k','linewidth',1);
% c = colorbar;
% myColorMap = jet; % Make a copy of jet.
% % Assign white (all 1's) to black (the first row in myColorMap).
% % myColorMap(1, :) = [1 1 1];
% colormap(myColorMap); % Apply the colormap


%   
% 
% load geoid  
% % Create a figure with an Eckert projection.  
% figure  
% axesm eckert4; %注意axesm后面的m了吗?,可以使用maps命令查看所有的地图投影的方式,然后选一个  
% framem; gridm;%显示框架和网格线,注意后面都多了个m,表示map  
% axis off %关闭外部坐标轴,外部坐标轴不同于map axes  
% % Display the geoid as a texture map.   
% geoshow(geoid, geoidrefvec, 'DisplayType', 'texturemap');  
% % Create a colorbar and title.  
% hcb = colorbar('southoutside');  
% set(get(hcb,'Xlabel'),'String','EGM96 Geoid Heights in Meters.')  
%  
% % Mask out all the land.  
% geoshow('landareas.shp', 'FaceColor', 'white');  
