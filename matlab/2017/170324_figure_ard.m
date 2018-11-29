figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
noOfLines = 16;
cmp = gray(noOfLines);

cmp(14,1)=0.83;
cmp(14,2)=0.83;
cmp(14,3)=0.83;
cmp(15,1)=0.86;
cmp(15,2)=0.86;
cmp(15,3)=0.86;
cmp(16,1)=0.9;
cmp(16,2)=0.9;
cmp(16,3)=0.9;

hLine = plot(x,plots);

for line = 1:noOfLines
set(hLine(line),'Color',cmp(line,:));
end


title('HWR=1.0, ес_{w}=0.4, е╓_{f}=0.3','FontSize',25,'fontweight','bold')

set(hLine,'LineWidth',3.5);
set(gca,'fontsize',18)

xlhand = get(gca,'xlabel')
set(xlhand,'string','Zenith Angle (deg.)','fontsize',22,'FontWeight','bold')

ylhand = get(gca,'ylabel')
set(ylhand,'string','Bulk Canyon Albedo','fontsize',22,'FontWeight','bold')

legend1 = legend('ес_{rd} = 0.05','0.06','0.07','0.08','0.09','0.10','0.11','0.12','0.13','0.14','0.15','0.16','0.17','0.18','0.19','0.20','Location','northwest')
set(legend1,'Position',[0.335763888888889 0.658004158004158 0.0828125 0.231115731115731]);
axis square