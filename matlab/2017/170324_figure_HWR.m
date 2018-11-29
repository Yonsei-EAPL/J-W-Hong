figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
noOfLines = 20;
cmp = gray(noOfLines);

cmp(17,1)=0.8;
cmp(17,2)=0.8;
cmp(17,3)=0.8;
cmp(18,1)=0.83;
cmp(18,2)=0.83;
cmp(18,3)=0.83;
cmp(19,1)=0.86;
cmp(19,2)=0.86;
cmp(19,3)=0.86;
cmp(20,1)=0.9;
cmp(20,2)=0.9;
cmp(20,3)=0.9;

hLine = plot(x,plots);

for line = 1:noOfLines
set(hLine(line),'Color',cmp(line,:));
end


title('ес_{rd}=0.1, ес_{w}=0.4, е╓_{f}=0.3','FontSize',25,'fontweight','bold')

set(hLine,'LineWidth',3.5);
set(gca,'fontsize',18)

xlhand = get(gca,'xlabel')
set(xlhand,'string','Zenith Angle (deg.)','fontsize',22,'FontWeight','bold')

ylhand = get(gca,'ylabel')
set(ylhand,'string','Bulk Canyon Albedo','fontsize',22,'FontWeight','bold')

legend1 = legend('HWR = 0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3','1.4','1.5','1.6','1.7','1.8','1.9','2.0','Location','northwest')
set(legend1,'Position',[0.335763888888889 0.658004158004158 0.0828125 0.231115731115731]);
axis square