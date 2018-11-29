x=zeros(1000,1);
y=x;
for i = 1:1000
    x(i,1) = 0.01*i;
%     y= 1- exp(-x/2);
y(i,1) = 0.001*exp(-x(i,1))+ 0.1*(1-exp(-x(i,1)));

end
clear i