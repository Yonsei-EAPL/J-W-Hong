%data_b = data(756001:792000,:);
%data_a = data(864001:900000,:);

%% 60_b
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data60_0714_b(:,5)); % date
v_bar = mean(data60_0714_b(:,6)); % date
w_bar = mean(data60_0714_b(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data60_0714_b(j,5)+sin(alpha)*data60_0714_b(j,6))+sin(beta)*data60_0714_b(j,7);
    temp_wind(j,2) = -sin(alpha)*data60_0714_b(j,5)+cos(alpha)*data60_0714_b(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data60_0714_b(j,5)+sin(alpha)*data60_0714_b(j,6))+cos(beta)*data60_0714_b(j,7);
    data60_0714_b(j,5) = temp_wind(j,1);
    data60_0714_b(j,6) = temp_wind(j,2);
    data60_0714_b(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data60_0714_b(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_60_0714_b = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data60_0714_b(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_60_0714_b = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data60_0714_b(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_60_0714_b = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data60_0714_b(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_60_0714_b = fpsd_s;
clear fpsd_s;

% C-spectrum
x=data60_0714_b(:,10);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_C_60_0714_b = fpsd_s;
clear fpsd_s;

% Q-spectrum
x=data60_0714_b(:,11);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_Q_60_0714_b = fpsd_s;
clear fpsd_s;
clear data60_0714_b

%% 60_a
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data60_0714_a(:,5)); % date
v_bar = mean(data60_0714_a(:,6)); % date
w_bar = mean(data60_0714_a(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data60_0714_a(j,5)+sin(alpha)*data60_0714_a(j,6))+sin(beta)*data60_0714_a(j,7);
    temp_wind(j,2) = -sin(alpha)*data60_0714_a(j,5)+cos(alpha)*data60_0714_a(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data60_0714_a(j,5)+sin(alpha)*data60_0714_a(j,6))+cos(beta)*data60_0714_a(j,7);
    data60_0714_a(j,5) = temp_wind(j,1);
    data60_0714_a(j,6) = temp_wind(j,2);
    data60_0714_a(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data60_0714_a(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_60_0714_a = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data60_0714_a(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_60_0714_a = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data60_0714_a(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_60_0714_a = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data60_0714_a(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_60_0714_a = fpsd_s;
clear fpsd_s;

% C-spectrum
x=data60_0714_a(:,10);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_C_60_0714_a = fpsd_s;
clear fpsd_s;

% Q-spectrum
x=data60_0714_a(:,11);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_Q_60_0714_a = fpsd_s;
clear fpsd_s;
clear data60_0714_a


%% 140_b
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data140_0714_b(:,5)); % date
v_bar = mean(data140_0714_b(:,6)); % date
w_bar = mean(data140_0714_b(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data140_0714_b(j,5)+sin(alpha)*data140_0714_b(j,6))+sin(beta)*data140_0714_b(j,7);
    temp_wind(j,2) = -sin(alpha)*data140_0714_b(j,5)+cos(alpha)*data140_0714_b(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data140_0714_b(j,5)+sin(alpha)*data140_0714_b(j,6))+cos(beta)*data140_0714_b(j,7);
    data140_0714_b(j,5) = temp_wind(j,1);
    data140_0714_b(j,6) = temp_wind(j,2);
    data140_0714_b(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data140_0714_b(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_140_0714_b = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data140_0714_b(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_140_0714_b = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data140_0714_b(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_140_0714_b = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data140_0714_b(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_140_0714_b = fpsd_s;
clear fpsd_s;
clear data140_0714_b


%% 140_a
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data140_0714_a(:,5)); % date
v_bar = mean(data140_0714_a(:,6)); % date
w_bar = mean(data140_0714_a(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data140_0714_a(j,5)+sin(alpha)*data140_0714_a(j,6))+sin(beta)*data140_0714_a(j,7);
    temp_wind(j,2) = -sin(alpha)*data140_0714_a(j,5)+cos(alpha)*data140_0714_a(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data140_0714_a(j,5)+sin(alpha)*data140_0714_a(j,6))+cos(beta)*data140_0714_a(j,7);
    data140_0714_a(j,5) = temp_wind(j,1);
    data140_0714_a(j,6) = temp_wind(j,2);
    data140_0714_a(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data140_0714_a(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_140_0714_a = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data140_0714_a(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_140_0714_a = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data140_0714_a(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_140_0714_a = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data140_0714_a(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_140_0714_a = fpsd_s;
clear fpsd_s;
clear data140_0714_a


%% 300_b
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data300_0714_b(:,5)); % date
v_bar = mean(data300_0714_b(:,6)); % date
w_bar = mean(data300_0714_b(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data300_0714_b(j,5)+sin(alpha)*data300_0714_b(j,6))+sin(beta)*data300_0714_b(j,7);
    temp_wind(j,2) = -sin(alpha)*data300_0714_b(j,5)+cos(alpha)*data300_0714_b(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data300_0714_b(j,5)+sin(alpha)*data300_0714_b(j,6))+cos(beta)*data300_0714_b(j,7);
    data300_0714_b(j,5) = temp_wind(j,1);
    data300_0714_b(j,6) = temp_wind(j,2);
    data300_0714_b(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data300_0714_b(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_300_0714_b = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data300_0714_b(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_300_0714_b = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data300_0714_b(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_300_0714_b = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data300_0714_b(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_300_0714_b = fpsd_s;
clear fpsd_s;

% C-spectrum
x=data300_0714_b(:,10);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_C_300_0714_b = fpsd_s;
clear fpsd_s;

% Q-spectrum
x=data300_0714_b(:,11);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_Q_300_0714_b = fpsd_s;
clear fpsd_s;
clear data300_0714_b


%% 300_a
% double rotation
temp_wind = zeros(36000,3);
u_bar = mean(data300_0714_a(:,5)); % date
v_bar = mean(data300_0714_a(:,6)); % date
w_bar = mean(data300_0714_a(:,7)); % date
alpha = atan2(v_bar,u_bar);
beta = atan2(w_bar,((u_bar^2 + v_bar^2)^(0.25)));
for j = 1:36000
    temp_wind(j,1) = cos(beta)*(cos(alpha)*data300_0714_a(j,5)+sin(alpha)*data300_0714_a(j,6))+sin(beta)*data300_0714_a(j,7);
    temp_wind(j,2) = -sin(alpha)*data300_0714_a(j,5)+cos(alpha)*data300_0714_a(j,6);
    temp_wind(j,3) = -sin(beta)*(cos(alpha)*data300_0714_a(j,5)+sin(alpha)*data300_0714_a(j,6))+cos(beta)*data300_0714_a(j,7);
    data300_0714_a(j,5) = temp_wind(j,1);
    data300_0714_a(j,6) = temp_wind(j,2);
    data300_0714_a(j,7) = temp_wind(j,3);
end
clear u_bar v_bar w_bar alpha beta j temp_wind

% u-spectrum
x=data300_0714_a(:,5);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_u_300_0714_a = fpsd_s;
clear fpsd_s;

% v-spectrum
x=data300_0714_a(:,6);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_v_300_0714_a = fpsd_s;
clear fpsd_s;

% w-spectrum
x=data300_0714_a(:,7);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_w_300_0714_a = fpsd_s;
clear fpsd_s;

% T-spectrum
x=data300_0714_a(:,8);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_T_300_0714_a = fpsd_s;
clear fpsd_s;

% C-spectrum
x=data300_0714_a(:,10);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_C_300_0714_a = fpsd_s;
clear fpsd_s;

% Q-spectrum
x=data300_0714_a(:,11);
nx=4;   % the number of variables
i=3;    %
dt=0.1;  % time difference between two conjecutive points (seconds)
nst=16; % starting point of non-overlapping average
pt_dec=10;  % how many points will be overlapped (point per decade)
mm=size(x);
n=mm(1);
n2=n/2;
nvar=mm(2);
xm(1,1)=mean(x);
xstd(1,1)=std(x);
xvar(1,1)=xstd(1,1).^2;
xp(:,1)=x(:,1)-xm(1,1); % perturbation value
fx(:,1)=fft(xp(:,1));
fxr(:,1)=real(fx(:,1))./real(n);
fxi(:,1)=imag(fx(:,1))./real(n);
s(:,1) = fxr(:,1).^2 + fxi(:,1).^2;
psd(:,1) = zeros(n2,1); % zero padding: power spectrum density
f = zeros(n2,1);  % zeros padding: natural frequency (1/sec)
for j=1:1:n2-1
    f(j) = real(j)/(real(n).*dt);
    psd(j,1) = (2.*n.*dt)*s(j+1,1)./xvar(1,1);  %
end
f(n2) = real(n2)/(real(n).*dt);
psd(n2,1) = (n.*dt)*s(n2+1,1)./xvar(1,1);
% display ('Test the Parseval Theorem')
% display ('The value below should be near around 1.0 for the energy conservation')
% % sum(psd)./(n.*dt)
fpsd(:,1)=zeros(n2,1); % je-woo 2012-11-17
fpsd(:,1) = f(:).*psd(:,1);
logx=zeros(n2,1); % je-woo 2012-11-17
for j = nst:1:n2
    logx(j) = log10(f(j));
end
lf_mi = logx(nst);
lf_ma = logx(n2);
nn = int8((lf_ma - lf_mi).*real(pt_dec));
ax=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    ax(js) = f(js);
end
fpsd_s=zeros(46,1); %(n2,1); % je-woo 2012-11-17
for js = 1:1:nst-1
    fpsd_s(js,1) = fpsd(js,1);
end
p = nst - 1;
k = nst - 1;
sp = 0.0;
for js = 1:1:45
    flag1 = 'F';
    lo = (lf_mi) + real(js-1)./real(pt_dec);
    hi = (lf_mi) + real(js)./real(pt_dec);
    fr = 0.0;
    index = 0;
    sp = 0.0;
    while (flag1 == 'F')
        p = p+1;
        if(logx(p)<hi)&&(p<n2)
            sp = sp + fpsd(p+1,1);
            fr = fr + f(p+1);
            index = index + 1;
            flah1 = 'F';
        else
            if(index>0)
                k = k+1;
                fpsd_s(k,1) = sp./real(index);
                ax(k) = fr./real(index);
            end
            p = p-1;
            flag1 = 'T';
        end
    end
end
clear ans dt f flag1 flah1 fpsd fr fx fxi fxr h hi i index j js k lf_ma
clear lf_mi lo logx mm n n2 nn nst nvar nx p psd pt_dec s sp xm xp xstd xvar %x
fpsd_s_Q_300_0714_a = fpsd_s;
clear fpsd_s;
clear data300_0714_a

%%
clear x

result0714 = zeros(46,32);
result0714(:,1) = fpsd_s_u_60_0714_b(:,1);
result0714(:,2) = fpsd_s_v_60_0714_b(:,1);
result0714(:,3) = fpsd_s_w_60_0714_b(:,1);
result0714(:,4) = fpsd_s_T_60_0714_b(:,1);
result0714(:,5) = fpsd_s_C_60_0714_b(:,1);
result0714(:,6) = fpsd_s_Q_60_0714_b(:,1);
clear fpsd_s_u_60_0714_b fpsd_s_v_60_0714_b fpsd_s_w_60_0714_b fpsd_s_T_60_0714_b fpsd_s_C_60_0714_b fpsd_s_Q_60_0714_b
result0714(:,7) = fpsd_s_u_60_0714_a(:,1);
result0714(:,8) = fpsd_s_v_60_0714_a(:,1);
result0714(:,9) = fpsd_s_w_60_0714_a(:,1);
result0714(:,10) = fpsd_s_T_60_0714_a(:,1);
result0714(:,11) = fpsd_s_C_60_0714_a(:,1);
result0714(:,12) = fpsd_s_Q_60_0714_a(:,1);
clear fpsd_s_u_60_0714_a fpsd_s_v_60_0714_a fpsd_s_w_60_0714_a fpsd_s_T_60_0714_a fpsd_s_C_60_0714_a fpsd_s_Q_60_0714_a
result0714(:,13) = fpsd_s_u_140_0714_b(:,1);
result0714(:,14) = fpsd_s_v_140_0714_b(:,1);
result0714(:,15) = fpsd_s_w_140_0714_b(:,1);
result0714(:,16) = fpsd_s_T_140_0714_b(:,1);
clear fpsd_s_u_140_0714_b fpsd_s_v_140_0714_b fpsd_s_w_140_0714_b fpsd_s_T_140_0714_b
result0714(:,17) = fpsd_s_u_140_0714_a(:,1);
result0714(:,18) = fpsd_s_v_140_0714_a(:,1);
result0714(:,19) = fpsd_s_w_140_0714_a(:,1);
result0714(:,20) = fpsd_s_T_140_0714_a(:,1);
clear fpsd_s_u_140_0714_a fpsd_s_v_140_0714_a fpsd_s_w_140_0714_a fpsd_s_T_140_0714_a
result0714(:,21) = fpsd_s_u_300_0714_b(:,1);
result0714(:,22) = fpsd_s_v_300_0714_b(:,1);
result0714(:,23) = fpsd_s_w_300_0714_b(:,1);
result0714(:,24) = fpsd_s_T_300_0714_b(:,1);
result0714(:,25) = fpsd_s_C_300_0714_b(:,1);
result0714(:,26) = fpsd_s_Q_300_0714_b(:,1);
clear fpsd_s_u_300_0714_b fpsd_s_v_300_0714_b fpsd_s_w_300_0714_b fpsd_s_T_300_0714_b fpsd_s_C_300_0714_b fpsd_s_Q_300_0714_b
result0714(:,27) = fpsd_s_u_300_0714_a(:,1);
result0714(:,28) = fpsd_s_v_300_0714_a(:,1);
result0714(:,29) = fpsd_s_w_300_0714_a(:,1);
result0714(:,30) = fpsd_s_T_300_0714_a(:,1);
result0714(:,31) = fpsd_s_C_300_0714_a(:,1);
result0714(:,32) = fpsd_s_Q_300_0714_a(:,1);
clear fpsd_s_u_300_0714_a fpsd_s_v_300_0714_a fpsd_s_w_300_0714_a fpsd_s_T_300_0714_a fpsd_s_C_300_0714_a fpsd_s_Q_300_0714_a
