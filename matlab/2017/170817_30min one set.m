%% model from 1Hz to 2Hz
 % because the range of bump is known as 1.39 Hz and around 1.39 Hz.

lx1 = 1.0044;
lx2 = 2.0039;
ly1 = 0.0377;
ly2 = 0.0279;

lx1 = log10(lx1);
lx2 = log10(lx2);
ly1 = log10(ly1);
ly2 = log10(ly2);

a = (ly2-ly1)/(lx2-lx1); %slope
b = ly2 - a*lx2; % offset

clear lx1 lx2 ly1 ly2

%%
tempx1 =0;
tempx2 =0;
for i = 1:9000
    if i>1
        if f(i,1)==1
            tempx1 = i
        elseif f(i,1)==2
            tempx2 = i-1
        end
    end
end
clear i

% nn= tempx2-tempx1+1;
% lx = zeros(nn,1);
% ly = zeros(nn,3);
% for i = tempx1:tempx2
%     lx(i-tempx1+1,1) = log10(f(i,1));
%     ly(i-tempx1+1,1) = log10(fpsd(i,1)); % log(original values)
%     ly(i-tempx1+1,2) = a*lx(i-tempx1+1,1) + b; % log(modelled values)
% end
% clear i
% fly = sum(ly(:,2))/sum(ly(:,1)); % factor for correction make original to modelled
% 
% for i = tempx1:tempx2
%     ly(i-tempx1+1,3) = fly*ly(i-tempx1+1,1); % corrected values
% end
% clear i
% sum(ly(:,2))/sum(ly(:,3))
% 
% 

fx2 = fx;
for i = 2:1799
    fx2(i,1)=0;
end
for i = 3600:14402
    fx2(i,1)=0;
end
for i = 16203:18000
    fx2(i,1)=0;
end
clear i
yy2 = ifft(fx2);


%%
for i = 1:18000
    data(i,10) = yy(i,1);
end
clear i 
for i = 1:18000
    data(i,10) = yy2(i,1);
end
clear i 

