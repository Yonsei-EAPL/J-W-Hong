figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','w');
noOfLines = 19;
cmp = gray(noOfLines);

cmp(16,1)=0.8;
cmp(16,2)=0.8;
cmp(16,3)=0.8;
cmp(17,1)=0.83;
cmp(17,2)=0.83;
cmp(17,3)=0.83;
cmp(18,1)=0.86;
cmp(18,2)=0.86;
cmp(18,3)=0.86;
cmp(19,1)=0.9;
cmp(19,2)=0.9;
cmp(19,3)=0.9;

hLine = plot(x,plots);

for line = 1:noOfLines
set(hLine(line),'Color',cmp(line,:));
end


title('HWR=1, ес_{rd}=0.1, ес_{w}=0.4','FontSize',25,'fontweight','bold')

set(hLine,'LineWidth',3.5);
set(gca,'fontsize',18)

xlhand = get(gca,'xlabel')
set(xlhand,'string','Zenith Angle (deg.)','fontsize',22,'FontWeight','bold')

ylhand = get(gca,'ylabel')
set(ylhand,'string','Bulk Canyon Albedo','fontsize',22,'FontWeight','bold')

legend1 = legend('е╓_{f} = 0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50','0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','Location','northwest')
set(legend1,'Position',[0.335763888888889 0.658004158004158 0.0828125 0.231115731115731]);
axis square