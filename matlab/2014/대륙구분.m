for i = 90:110
    for j = 60:80
        if land(i,j)>0
            land(i,j) = 2; % east asia
        end
    end
end
clear i j

for i = 20:85
    for j = 140:180
        if land(i,j)>0
            land(i,j) = 3; % south america
        end
    end
end
clear i j

for i = 40:104
    for j = 1:28
        if land(i,j)>0
            land(i,j) = 4; % africa
        end
    end
end
clear i j
for i = 60:111
    for j = 180:192
        if land(i,j)>0
            land(i,j) = 4; % africa
        end
    end
end
clear i j

for i = 65:89
    for j = 50:90
        if land(i,j)>0
            land(i,j) = 5; % south-asia
        end
    end
end
clear i j

for i = 30:64
    for j = 60:100
        if land(i,j)>0
            land(i,j) = 6; % austrailia
        end
    end
end
clear i j

for i = 86:145
    for j = 103:176
        if land(i,j)>0
            land(i,j) = 7; % north-america
        end
    end
end
clear i j

for i = 105:145
    for j = 1:40
        if land(i,j)>0
            land(i,j) = 8; % europe
        end
    end
end
clear i j
for i = 112:145
    for j = 178:192
        if land(i,j)>0
            land(i,j) = 8; % europe
        end
    end
end
clear i j

for i = 1:145
    for j = 1:192
        if (land(i,j)>0)&&(land(i,j)<2)
            land(i,j) = 1;
        end
    end
end
clear i j