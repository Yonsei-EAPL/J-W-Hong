%%
fra1 = zeros(1001,1001);
fra2 = zeros(1001,1001);
fra3 = zeros(1001,1001);

for i = 1:1001
    for j = 1:1001
        fra1(i,j) = fraction3(i,j,1);
        fra2(i,j) = fraction3(i,j,2);
        fra3(i,j) = fraction3(i,j,3);
    end
end
clear i j

for i = 1:1001
    for j = 1:1001
        if ((fra1(i,j) ==0)&&(fra2(i,j)==0))&&(fra3(i,j)==0)
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==224))&&(fra3(i,j)==0)
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==192))&&(fra3(i,j)==192)
        elseif ((fra1(i,j) ==128)&&(fra2(i,j)==128))&&(fra3(i,j)==128)
        elseif ((fra1(i,j) ==195)&&(fra2(i,j)==195))&&(fra3(i,j)==195)
            fra1(i,j) =192;
            fra2(i,j) =192;
            fra3(i,j) =192;            
        else
            fra1(i,j) =192;
            fra2(i,j) =224;
            fra3(i,j) =0;
        end
    end
end
clear i j

fra = zeros(1001,1001,3);
for i = 1:1001
    for j = 1:1001
        fra(i,j,1) = fra1(i,j)/255;
        fra(i,j,2) = fra2(i,j)/255;
        fra(i,j,3) = fra3(i,j)/255;
    end
end
clear i j

imshow(fra)


%%
lulc = zeros(1001,1001);

for i = 1:1001
    for j = 1:1001
        if ((fra1(i,j) ==192)&&(fra2(i,j)==224))&&(fra3(i,j)==0)
            lulc(i,j) = 1; % veg
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==192))&&(fra3(i,j)==192)
            lulc(i,j) = 2; % building
        else
            lulc(i,j) = 3; % road
        end
        
    end
end
clear i j 

lulc_wd = zeros(72,4); % 5-degree interval 

x=501;
y=501;
for i = 1:1001
    for j = 1:1001
        if ((i-x)^2+(j-y)^2)^0.5 <=500
            ang = atan2(i-x,j-y)/pi()*180;
            if ang<0
                ang=ang+360;
            end
            ang = (ang-mod(ang,5))/5;
            if ang ==0
                ang =72;
            end
            if lulc(i,j)==1
                lulc_wd(ang,1) = lulc_wd(ang,1)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;
            elseif lulc(i,j)==2
                lulc_wd(ang,2) = lulc_wd(ang,2)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;            
            else
                lulc_wd(ang,3) = lulc_wd(ang,3)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;            
            end
        end
    end
end
clear i j

lulc_wd(73,1)= sum(lulc_wd(1:72,1));
lulc_wd(73,2)= sum(lulc_wd(1:72,2));
lulc_wd(73,3)= sum(lulc_wd(1:72,3));
lulc_wd(73,4)= sum(lulc_wd(1:72,4));

lulc_wd_fra = lulc_wd;
for i = 1:73
    for j = 1:4
        lulc_wd_fra(i,j) = lulc_wd_fra(i,j)/lulc_wd_fra(i,4)*100;
    end
end
clear i j












%%
fra1 = zeros(1001,1001);
fra2 = zeros(1001,1001);
fra3 = zeros(1001,1001);

for i = 1:1001
    for j = 1:1001
        fra1(i,j) = fraction4(i,j,1);
        fra2(i,j) = fraction4(i,j,2);
        fra3(i,j) = fraction4(i,j,3);
    end
end
clear i j

for i = 1:1001
    for j = 1:1001
        if ((fra1(i,j) ==0)&&(fra2(i,j)==0))&&(fra3(i,j)==0)
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==224))&&(fra3(i,j)==0)
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==192))&&(fra3(i,j)==192)
        elseif ((fra1(i,j) ==0)&&(fra2(i,j)==162))&&(fra3(i,j)==232)
        elseif ((fra1(i,j) ==128)&&(fra2(i,j)==128))&&(fra3(i,j)==128)
        elseif ((fra1(i,j) ==195)&&(fra2(i,j)==195))&&(fra3(i,j)==195)
            fra1(i,j) =192;
            fra2(i,j) =192;
            fra3(i,j) =192;            
        else
            fra1(i,j) =192;
            fra2(i,j) =224;
            fra3(i,j) =0;
        end
    end
end
clear i j

fra = zeros(1001,1001,3);
for i = 1:1001
    for j = 1:1001
        fra(i,j,1) = fra1(i,j)/255;
        fra(i,j,2) = fra2(i,j)/255;
        fra(i,j,3) = fra3(i,j)/255;
    end
end
clear i j

imshow(fra)

%%
lulc = zeros(1001,1001);

for i = 1:1001
    for j = 1:1001
        if ((fra1(i,j) ==192)&&(fra2(i,j)==224))&&(fra3(i,j)==0)
            lulc(i,j) = 1; % veg
        elseif ((fra1(i,j) ==192)&&(fra2(i,j)==192))&&(fra3(i,j)==192)
            lulc(i,j) = 2; % building
        elseif ((fra1(i,j) ==0)&&(fra2(i,j)==162))&&(fra3(i,j)==232)
            lulc(i,j) = 4; % research-building (2 times of building for pavement sfc.)
        else
            lulc(i,j) = 3; % road
        end
        
    end
end
clear i j 

lulc_wd = zeros(72,4); % 5-degree interval 

x=501;
y=501;
for i = 1:1001
    for j = 1:1001
%         if ((i-x)^2+(j-y)^2)^0.5 <=500
            ang = atan2(i-x,j-y)/pi()*180;
            if ang<0
                ang=ang+360;
            end
            ang = (ang-mod(ang,5))/5;
            if ang ==0
                ang =72;
            end
            if lulc(i,j)==1
                lulc_wd(ang,1) = lulc_wd(ang,1)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;
            elseif lulc(i,j)==2
                lulc_wd(ang,2) = lulc_wd(ang,2)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;            
            elseif lulc(i,j)==4
                lulc_wd(ang,2) = lulc_wd(ang,2)+1;
                lulc_wd(ang,3) = lulc_wd(ang,3)+1;
                lulc_wd(ang,1) = lulc_wd(ang,1)-1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;
            else
                lulc_wd(ang,3) = lulc_wd(ang,3)+1;
                lulc_wd(ang,4) = lulc_wd(ang,4)+1;            
            end
%         end
    end
end
clear i j

lulc_wd(73,1)= sum(lulc_wd(1:72,1));
lulc_wd(73,2)= sum(lulc_wd(1:72,2));
lulc_wd(73,3)= sum(lulc_wd(1:72,3));
lulc_wd(73,4)= sum(lulc_wd(1:72,4));

lulc_wd_fra = lulc_wd;
for i = 1:73
    for j = 1:4
        lulc_wd_fra(i,j) = lulc_wd_fra(i,j)/lulc_wd_fra(i,4)*100;
    end
end
clear i j
