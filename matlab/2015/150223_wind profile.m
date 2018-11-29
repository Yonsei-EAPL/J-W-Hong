for i = 1:144
    hold on
    plot(data(i,:),y(1,:),'x-b','LineWidth',4,'MarkerSize',20)
    hold on
    plot(data2(i,:),y2(1,:),'k')
end
clear i
