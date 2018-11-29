%%
% 1<- 102 ; building 
% 2<- 0 ; road
% 3<- 164 ; impervious
% 4<- 210 ; water
% 5<- 113 ; forest
% 6<- 62 ; grass
% 7<- 8 ; garden
% 8<- 251 ; bare-soil


%%
LULC = cdata;
clear cdata

for i = 1:999
    for j = 1:999
        if LULC(i,j)==102
            LULC(i,j) = 1;
        elseif LULC(i,j)==0
            LULC(i,j) = 2;
        elseif LULC(i,j)==164
            LULC(i,j) = 3;
        elseif LULC(i,j)==210
            LULC(i,j) = 4;
        elseif LULC(i,j)==113
            LULC(i,j) = 5;
        elseif LULC(i,j)==62
            LULC(i,j) = 6;
        elseif LULC(i,j)==8
            LULC(i,j) = 7;
        elseif LULC(i,j)==251
            LULC(i,j) = 8;
        end
    end
end
clear i j

