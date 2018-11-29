%% building cover
%   source : 'H:\b1_문서\b1_doc_2012_대학원\1302_은평뉴타운\140702_land cover\lulc3_bmp.bmp'
%   simply devide landcover as building and road only
bd = cdata(501:999,:);
clear cdata colormap
for i = 1:499
    for j = 1:999
        if bd(i,j)==102
            bd(i,j)=20; % building
        else
            bd(i,j)=0;  % road
        end
    end
end
clear i j

%% 90 deg cover
%   only use 135-225 deg for footprint analysis area of radiometer(CNR4)
%   and within 500-m distance from tower (0, 500) << which is not included
%   in bd & bc
dc = zeros(499,999);
for i = 1:499
    for j = 1:999
        if (i^2 + (j-500)^2)^0.5<500
            if abs(atan((500-j)/i))<(atan(1))
                dc(i,j)=1;
            end
        end
    end
end
clear i j 

%% make footprint
%   F - col1 : accumulated footprint(%), col2 : distance (r, m)
F(:,3)=0;
fc1 = zeros(499,999);
for i = 1:499
    for j = 1:999
        if dc(i,j)==1
            r = ((i^2 + (j-500)^2)^0.5);
            if r <= F(19,2);
                temp1=0;
                temp2=0;
                for k = 1:19
                    if (F(k,2)>r)&&(temp1==0)
                        temp1 = 1;
                        temp2 = k;
                    end
                end
                clear k
                F(temp2,3) = F(temp2,3)+1;
                fc1(i,j) = temp2;
            else
                F(20,3) = F(20,3)+1;
                fc1(i,j)=20;
            end
        end
    end
end
clear i j r

F(:,4)=0;
for i = 1:19
    F(i,4) = 5/F(i,3);
end
clear i

%% footprint cover
fc2 = zeros(499,999);
for i = 1:499
    for j = 1:999
        if fc1(i,j)>=1;
            fc2(i,j) = F(fc1(i,j),4);
        end
    end
end
clear i j

%% building fraction
result_BF = 0;
for i = 1:499
    for j = 1:999
        if (bd(i,j)>0)&&(fc2(i,j)>0)
            result_BF = result_BF+fc2(i,j);
        end
    end
end
clear i j
