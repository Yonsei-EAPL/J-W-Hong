x=0:0.01:24*pi;
y1=sin(x);
y2=sin(x+pi/2);
[r,lags] = xcov(y1(:),y2(:));
figure
plot(lags,r,'-o');
[Y,I] = max(r);
Y, I=I-length(x)
[Y,I] = min(r);
Y, I=I-length(x)
