pano = zeros(500,1440,3); % 500m * 0.25 degree * 3 colors
angle = zeros(1039,1120);
angle2 = zeros(1039,1120);
angle3 = zeros(1039,1120);
temp = 0;
for i = 21:1020
    for j = 71:1070

        temp = atan2(i-520,j-570)/pi()*180;
        temp = temp+90; % last
        if temp<0
            temp=temp+360;
        end
        if temp == 0
            temp = 360;
        end
        angle(i,j) = temp;
        x = (temp - mod(temp,0.25))/0.25;
        if x>0
            y = ((j-570)^2+(i-520)^2)^(0.5);
            y = round(y);
            angle2(i,j) = y;
            angle3(i,j) = x;
            if (y<=500)&&(y>0)
                if pano(501-y,x,1)==0
                    for z = 1:3
                        pano(501-y,x,z) = landscape2(i,j,z);
                    end
                    clear z
                end
            end
        end  
    end
end
clear i j

for i = 1:500
    for j = 1:1440
        for z = 1:3
            pano(i,j,z) =pano(i,j,z)/255;
        end
    end
end
clear i j z

imshow(pano)


%% interpolation

temp=0;
for i = 1:500
    for j = 1:1440
        temp2=0;
        for k = 1:3
            temp2 = temp2 + pano(i,j,k);
        end
        if temp2~=0
            temp = temp+1;
        end
    end
end
clear i j k temp2

panoxyz = zeros(temp,5);
temp=0;
for i = 1:500
    for j = 1:1440
        temp2=0;
        for k = 1:3
            temp2 = temp2 + pano(i,j,k);
        end
        if temp2~=0
            temp = temp+1;
            panoxyz(temp,1) = i;
            panoxyz(temp,2) = j;
            panoxyz(temp,3) = pano(i,j,1);
            panoxyz(temp,4) = pano(i,j,2);
            panoxyz(temp,5) = pano(i,j,3);
        end
    end
end
clear i j k temp2


X = 1:1:1440;
Y = 1:1:500;
[XI, YI] = meshgrid(X,Y);
ZI1 = griddata(panoxyz(:,2),panoxyz(:,1),panoxyz(:,3),XI,YI);
ZI2 = griddata(panoxyz(:,2),panoxyz(:,1),panoxyz(:,4),XI,YI);
ZI3 = griddata(panoxyz(:,2),panoxyz(:,1),panoxyz(:,5),XI,YI);

panoxyz= zeros(500,1440,3);
for i = 1:500
    for j = 1:1440
        panoxyz(i,j,1) = ZI1(i,j);
        panoxyz(i,j,2) = ZI2(i,j);
        panoxyz(i,j,3) = ZI3(i,j);
    end
end
clear i j



%%

% pano = zeros(500,360,3); % 500m * 1 degree * 3 colors
% angle = zeros(1039,1120);
% angle2 = zeros(1039,1120);
% angle3 = zeros(1039,1120);
% temp = 0;
% for i = 21:1020
%     for j = 71:1070
% 
%         temp = atan2(i-520,j-570)/pi()*180;
%         if temp<0
%             temp=temp+360;
%         end
%         if temp == 0
%             temp = 360;
%         end
%         angle(i,j) = temp;
%         x = round(temp);
%         if x>0
%             y = ((j-570)^2+(i-520)^2)^(0.5);
%             y = round(y);
%             angle2(i,j) = y;
%             angle3(i,j) = x;
%             if (y<=500)&&(y>0)
%                 if pano(501-y,x,1)==0
%                     for z = 1:3
%                         pano(501-y,x,z) = landscape2(i,j,z);
%                     end
%                     clear z
%                 end
%             end
%         end  
%     end
% end
% clear i j
% 
% for i = 1:500
%     for j = 1:360
%         for z = 1:3
%             pano(i,j,z) =pano(i,j,z)/255;
%         end
%     end
% end
% clear i j z
% 
% imshow(pano)
% 
