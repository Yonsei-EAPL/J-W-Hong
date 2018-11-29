%% wavelet
    % temp(:,9) - horizontal windspeed, 
    % temp(:,10) - total windspeed, 

% for i = 1:18000
%     temp(i,9) = (temp(i,1)^2 + temp(i,2)^2)^(0.5);
%     temp(i,10) = (temp(i,1)^2 + temp(i,2)^2 + temp(i,3)^2)^(0.5);
% end
% clear i 


%% figure
    % time-series
% figure()
% plot(temp(:,1))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('u','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('u (m s-1)','fontweight','bold','FontSize',20);
% figure()
% plot(temp(:,2))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('v','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('v (m s-1)','fontweight','bold','FontSize',20);
% figure()
% plot(temp(:,9))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('horizontal windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('horizontal windspeed (m s-1)','fontweight','bold','FontSize',20);
% figure()
% plot(temp(:,10))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('total windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('windspeed (m s-1)','fontweight','bold','FontSize',20);
% figure()
% plot(temp(:,7))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('amb-pressure','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('amb-pressure (kPa)','fontweight','bold','FontSize',20);    
%     



    % wt
% figure()
% wt(temp(:,1))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('WT of u','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wt(temp(:,2))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('WT of v','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wt(temp(:,9))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('WT of horizontal windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wt(temp(:,10))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('WT of total windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wt(temp(:,7))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('WT of amb-pressure','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);

    % wtc (coherence)
% figure()
% wtc(temp(:,1),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of u','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wtc(temp(:,2),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&v','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wtc(temp(:,9),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&horizontal windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wtc(temp(:,10),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&total windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% wtc(temp(:,7),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&amb-pressure','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);

    % xwt (cross)
% figure()
% xwt(temp(:,1),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of u','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% xwt(temp(:,2),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&v','fontweight','bold','FontSize',28)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% xwt(temp(:,9),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&horizontal windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% xwt(temp(:,10),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&total windspeed','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
% figure()
% xwt(temp(:,7),temp(:,5))
% set(gcf,'color','w');
% set(gca,'fontweight','bold','FontSize',20);
% title('Cross-WT of co2&amb-pressure','fontweight','bold','FontSize',24)
% xlabel('n (10Hz)','fontweight','bold','FontSize',20);
% ylabel('Period (n)','fontweight','bold','FontSize',20);
