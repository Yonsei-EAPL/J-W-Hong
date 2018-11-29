figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
noOfLines = 7;
cmp = gray(noOfLines);

cmp(7,1)=0.9;
cmp(7,2)=0.9;
cmp(7,3)=0.9;

hLine = plot(x,plots);

for line = 1:noOfLines
set(hLine(line),'Color',cmp(line,:));
end


title('HWR=1.0, ес_{rd}=0.1, е╓_{f}=0.3','FontSize',25,'fontweight','bold')

set(hLine,'LineWidth',3.5);
set(gca,'fontsize',18)

xlhand = get(gca,'xlabel')
set(xlhand,'string','Zenith Angle (deg.)','fontsize',22,'FontWeight','bold')

ylhand = get(gca,'ylabel')
set(ylhand,'string','Bulk Canyon Albedo','fontsize',22,'FontWeight','bold')

legend1 = legend('ес_w = 0.20','0.25','0.30','0.35','0.40','0.45','0.50','Location','northwest')
set(legend1,'Position',[0.335763888888889 0.658004158004158 0.0828125 0.231115731115731]);
axis square