%% total analysis
result = zeros(8,2);    % total
result(1,1) = 0;    % roads
result(2,1) = 8;    % garden
result(3,1) = 62;   % grass
result(4,1) = 102;  % buildings
result(5,1) = 113;  % forest
result(6,1) = 164;  % impervious(pavements)
result(7,1) = 210;  % water
result(8,1) = 251;  % baresoil
a=0;    % error index
for i = 301:700 %1:999
    for j = 301:700 %1:999
        if landcover(i,j)==0
            result(1,2)=result(1,2)+1;
        elseif landcover(i,j)==8
            result(2,2)=result(2,2)+1;
        elseif landcover(i,j)==62
            result(3,2)=result(3,2)+1;
        elseif landcover(i,j)==102
            result(4,2)=result(4,2)+1;
        elseif landcover(i,j)==113
            result(5,2)=result(5,2)+1;
        elseif landcover(i,j)==164
            result(6,2)=result(6,2)+1;
        elseif landcover(i,j)==210
            result(7,2)=result(7,2)+1;
        elseif landcover(i,j)==251
            result(8,2)=result(8,2)+1;
        else
            a=a+1
        end
    end
end
clear i j a


%% directional analysis
for i = 301:700 %1:999
    for j = 301:700 %1:999
        if direction(i,j) == 79
            direction(i,j) = 2;
        elseif direction(i,j) == 210
            direction(i,j) = 3;
        elseif direction(i,j) == 103
            direction(i,j) = 4;            
        elseif direction(i,j) == 232
            direction(i,j) = 5;            
        elseif direction(i,j) == 251
            direction(i,j) = 6;
        elseif direction(i,j) == 113
            direction(i,j) = 7;
        elseif direction(i,j) == 164
            direction(i,j) = 8;
        elseif direction(i,j) == 255
            direction(i,j) = 9;
        end
    end
end
clear i j 
result_dir = zeros(8,9);    %  8-direction (from NNE'0-45' to NNW'315-360')
result_dir(1,1) = 0;    % roads
result_dir(2,1) = 8;    % garden
result_dir(3,1) = 62;   % water
result_dir(4,1) = 102;  % buildings
result_dir(5,1) = 113;  % forest
result_dir(6,1) = 164;  % impervious(pavements)
result_dir(7,1) = 210;  % grass
result_dir(8,1) = 251;  % baresoil
for i = 301:700 %1:999
    for j = 301:700 %1:999
        if landcover(i,j)==0
            result_dir(1,direction(i,j))=result_dir(1,direction(i,j))+1;
        elseif landcover(i,j)==8
            result_dir(2,direction(i,j))=result_dir(2,direction(i,j))+1;
        elseif landcover(i,j)==62
            result_dir(3,direction(i,j))=result_dir(3,direction(i,j))+1;
        elseif landcover(i,j)==102
            result_dir(4,direction(i,j))=result_dir(4,direction(i,j))+1;
        elseif landcover(i,j)==113
            result_dir(5,direction(i,j))=result_dir(5,direction(i,j))+1;
        elseif landcover(i,j)==164
            result_dir(6,direction(i,j))=result_dir(6,direction(i,j))+1;
        elseif landcover(i,j)==210
            result_dir(7,direction(i,j))=result_dir(7,direction(i,j))+1;
        elseif landcover(i,j)==251
            result_dir(8,direction(i,j))=result_dir(8,direction(i,j))+1;
        end
    end
end
clear i j





