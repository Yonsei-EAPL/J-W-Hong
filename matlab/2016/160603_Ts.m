y=fft(x);
yreal = real(y);
yimag = imag(y);

for i = 1:2889
    yimag(i,1) = yimag(i,1)^2;
    yreal(i,1) = yreal(i,1)^2;
end
clear i

Sy = zeros(1444,1);
for i = 1:1444
    Sy(i,1) = yimag(i+1,1)+yimag(2890-i,1)+yreal(i+1,1)+yreal(2890-i,1);
end
clear i

time = double(time)/24;
