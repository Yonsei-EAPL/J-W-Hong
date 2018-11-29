
figure()
loglog(ax(:),result0714(:,6),'-ob')
hold on
loglog(ax(:),result0730(:,6),'-ob')
loglog(ax(:),result0801(:,6),'-ob')
loglog(ax(:),result0803(:,6),'-ob')
loglog(ax(:),result0806(:,6),'-ob')
loglog(ax(:),result0815(:,6),'-ob')
loglog(ax(:),result0908(:,6),'-ob')
loglog(ax(:),result0909(:,6),'-ob')
loglog(ax(:),result0910(:,6),'-ob')

loglog(ax(:),result0714(:,12),'-or')
hold on
loglog(ax(:),result0730(:,12),'-or')
loglog(ax(:),result0801(:,12),'-or')
loglog(ax(:),result0803(:,12),'-or')
loglog(ax(:),result0806(:,12),'-or')
loglog(ax(:),result0815(:,12),'-or')
loglog(ax(:),result0908(:,12),'-or')
loglog(ax(:),result0909(:,12),'-or')
loglog(ax(:),result0910(:,12),'-or')





x = data(468001:1296000,8);
wt(x)