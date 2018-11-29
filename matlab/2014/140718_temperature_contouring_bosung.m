%% History
% 140717 Mr.Keunmin Lee; code

%% Make contour field
% First, make array for ploting contour
% x(time:1~48), y(height: 10,20,40,60,80,100,140,180,220,260,300), z(mean_Temperature)
% input : mean_temp

%mean_WS=sqrt(mean_WD_u.^2+mean_WD_v.^2);
height = [10,20,40,60,80,100,140,180,220,260,300];
k=0;
for i=1:48
    for j=1:11
        k=k+1;
        x(k)=i;
        y(k)=height(j);        
        z(k)=mean_temp(i,j);
    end
end
clear i j k height

% Second, define a regular grid. I'll set up a 291x48 grid. [XI,YI]
x_step=1:1:48;
y_step=1:1:300;
[XI,YI] = meshgrid(x_step, y_step);
% XI contains the x-coords of each point and YI contains the y-coords.
ZI = griddata(x,y,z,XI,YI);
clear x y z x_step y_step XI YI

% plot 'imagesc'
imagesc(ZI);figure(gcf);
set(gca,'YDir','normal');
axis([1 48 10 300])
set(gcf,'color','w');
title('JUNE mean diurnal variation of Temperature')
xlabel('Time (30min mean)')
ylabel('Height(m)')
set(gca,'XTickLabel',['02:30';'05:00';'07:30';'10:00';'12:30';'15:00';'17:30';'20:00';'22:30'])

