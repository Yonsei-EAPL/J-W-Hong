t = 1 : 1000; t = t';
x = 2*sin(2*pi*t/50) + sin(2*pi*t/15) +0.5*sin(2*pi*t/5);
randn('seed',0)
n = randn(1000,1);
xn = x + n;
xt = x + 0.005*t;
Xxx = fft(x,1024);


%%
Pxx2 = abs(Xxx).^2/1000;
Pxx = [Pxx2(1); 2*Pxx2(2:512)];
f = 0 : 1/(1024-1) : 1/2;
subplot(1,2,1)
plot(f,Pxx), grid
axis([0 0.5 0 max(Pxx)])
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')
%%
[Pxx,f] = periodogram(x,[],1024,1);
subplot(1,2,2)
plot(f,Pxx), grid
axis([0 0.5 0 max(Pxx)])
xlabel('Frequency')
ylabel('Power')
title('Auto-Spectrum')

clear
t = 1 : 1000; t = t';
x = 2*sin(2*pi*t/50)+sin(2*pi*t/15)+0.5*sin(2*pi*t/5);
randn('seed',0)
n = randn(1000,1);
xn = x + n;
xt = x + 0.005*t;
xdt = detrend(xt);
subplot(2,1,1)
plot(t,x,'b-',t,xt,'r-'), grid
axis([0 200 -4 4])
subplot(2,1,2)
plot(t,x,'b-',t,xdt,'r-'), grid
axis([0 200 -4 4])
[Pxxt,f] = periodogram(xt,[],1024,1);
[Pxxdt,f]=periodogram(xdt,[],1024,1);
figure()
subplot(1,2,1)
plot(f,Pxxt), grid
xlabel('Frequency')
ylabel('Power')
subplot(1,2,2)
plot(f,Pxxdt), grid
xlabel('Frequency')
ylabel('Power')

clear
t = 1 : 1000;
% x = 2*sin(2*pi*t/5);
y = 2*sin(2*pi*t/5 + 2*pi/5);
x = 2*sin(2*pi*t/10);
% y = 2*sin(2*pi*t/10 + 2*pi/5);
figure()
plot(t,x,'b-',t,y,'r-')
axis([0 50 -2 2]), grid
[Pxy,f] = cpsd(x,y,[],0,1024,1);
figure()
plot(f,abs(Pxy)), grid
xlabel('Frequency')
ylabel('Power')
title('Cross-Spectrum')


clear
t = 1 : 1000;
x = 2*sin(2*pi*t/5);
%y = 2*sin(2*pi*t/5 +2*pi/5)+0.1*randn(size(t));
% y = 2*sin(2*pi*t/5 +1*pi/2)+0.1*randn(size(t));
y = 2*sin(2*pi*t/5 -2*pi/5);%+0.1*randn(size(t));
%y = 2*sin(2*pi*t/5 );%+0.1*randn(size(t));
figure()
plot(t,x,'b-',t,y,'r-')
axis([0 50 -2 2]), grid
[Pxy,f] = cpsd(x,y,[],0,1024,1);
figure()
subplot(2,1,1)
plot(f,abs(Pxy)), grid
% xlabel('Frequency')
% ylabel('Power')
% title('Cross-Spectrum')
% [Cxy,f] =
% %mscohere(x,y,[],0,1024,1);
% plot(f,Cxy), grid
% xlabel('Frequency')
% ylabel('Coherence')
% title('Coherence')
subplot(2,1,2)
phase = angle(Pxy);
plot(f,phase), grid
xlabel('Frequency')
ylabel('Phase Angle')
title('Phase Spectrum')




%% FFT
clear
t=1:1000; t=t';
x=2*sin(2*pi*t/50)+sin(2*pi*t/15)+0.5*sin(2*pi*t/5);
randn('seed',0)
n=randn(1000,1);
xn = x+n;
xt = x+0.005*t;

Xxx = fft(x,1024);
Pxx2 = abs(Xxx).^2/1000;
Pxx = [Pxx2(1); 2*Pxx2(2:512)];
f=0:1/(1024-1):1/2;
figure()
subplot(1,2,1)
plot(f,Pxx),grid
axis([0 0.5 0 max(Pxx)])
xlabel('Frequency')
ylabel('Power')
title('Auto-Apectrum')

%% Periodogram
[Period_Pxx, Period_f] = periodogram(x,[],1024,1);
subplot(1,2,2)
plot(Period_f,Period_Pxx),grid
axis([0 0.5 0 max(Period_Pxx)])

xlabel('Frequency')
ylabel('Power')
title('Auto-Apectrum')




