%% histogram for result

figure(1)

subplot(2,2,1)
hist(result2006(:,2),50)
title('AWS 2006')
hold on

subplot(2,2,2)
hist(result2007(:,2),50)
title('AWS 2007')
hold on

subplot(2,2,3)
hist(result2008(:,2),50)
title('AWS 2008')
hold on

subplot(2,2,4)
hist(result2012(:,2),50)
title('AWS 2012')
hold on

