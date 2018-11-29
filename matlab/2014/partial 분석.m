%%

% x1 = zeros(360*10,145,192); %isoprene; daily mean (BVOC.nc)
% x2 = zeros(360*10,145,192); %gpp; daily mean (GPP.nc)
% x3 = zeros(360*10,145,192); %sw_down; daily mean (Gu.nc)
% x4 = zeros(360*10,145,192); %q_air; daily mean (drive)
% x5 = zeros(360*10,145,192); %t_air; daily mean (drive)
% x6 = zeros(360*10,145,192); % precip; daily accumulated (drive)


%%

% yr = 10;
% 
% iso = nc_varget('control.BVOC.nc','isoprene_gb');
% gpp = nc_varget('control.GPP.nc','gpp_gb');
% sw = nc_varget('control.Gu.nc','sw_down');
% 
% for i = 1:360
%     for j = 1:8
%         for m = 1:7246
%             x1((yr-1)*360+i,add(1,m),add(2,m)) = x1((yr-1)*360+i,add(1,m),add(2,m))+iso((i-1)*8+j,m)*10800*1000;
%             x2((yr-1)*360+i,add(1,m),add(2,m)) = x2((yr-1)*360+i,add(1,m),add(2,m))+gpp((i-1)*8+j,m)*10800*1000;
%             x3((yr-1)*360+i,add(1,m),add(2,m)) = x3((yr-1)*360+i,add(1,m),add(2,m))+sw((i-1)*8+j,m)*10800/(10^6);
%         end
%     end
% end
% clear i j m
% clear iso gpp sw


%%

% yr = 10;
% 
% q = nc_varget('Qair_AMIP_2005.nc','Qair');
% t = nc_varget('Tair_AMIP_2005.nc','Tair');
% pr = nc_varget('Rainf_AMIP_2005.nc','Rainf');
% ps = nc_varget('Snowf_AMIP_2005.nc','Snowf');
% p = pr+ps;
% clear pr ps
% 
% for i = 1:360
%     for j = 1:8
%         for m = 1:145
%             for n = 1:192
%                 x4((yr-1)*360+i,m,n) = x4((yr-1)*360+i,m,n) + q((i-1)*8+j,m,n)/8;
%                 x5((yr-1)*360+i,m,n) = x5((yr-1)*360+i,m,n) + t((i-1)*8+j,m,n)/8;
%                 x6((yr-1)*360+i,m,n) = x6((yr-1)*360+i,m,n) + p((i-1)*8+j,m,n)*10800;
%             end
%         end
%     end
% end
% clear i j m n
% clear q t p


%%
clear yr
x1(~isfinite(x1))=0;
x2(~isfinite(x2))=0;
x3(~isfinite(x3))=0;
x4(~isfinite(x4))=0.05;
x5(~isfinite(x5))=273;
x6(~isfinite(x6))=0;

%%
pc12 = zeros(145,192);
pc13 = zeros(145,192);
pc14 = zeros(145,192);
pc15 = zeros(145,192);
pc16 = zeros(145,192);
for i = 1:145
    for j = 1:192
        if land(i,j)>0
            pc12(i,j) = partialcorr6(x1(:,i,j),x2(:,i,j),x3(:,i,j),x4(:,i,j),x5(:,i,j),x6(:,i,j));
            pc13(i,j) = partialcorr6(x1(:,i,j),x3(:,i,j),x2(:,i,j),x4(:,i,j),x5(:,i,j),x6(:,i,j));
            pc14(i,j) = partialcorr6(x1(:,i,j),x4(:,i,j),x2(:,i,j),x3(:,i,j),x5(:,i,j),x6(:,i,j));
            pc15(i,j) = partialcorr6(x1(:,i,j),x5(:,i,j),x2(:,i,j),x3(:,i,j),x4(:,i,j),x6(:,i,j));
            pc16(i,j) = partialcorr6(x1(:,i,j),x6(:,i,j),x2(:,i,j),x3(:,i,j),x4(:,i,j),x5(:,i,j));
        end
    end
end
clear i j

for i = 1:145
    for j = 1:192
        if (pc12(i,j)>-0.2)&&(pc12(i,j)<0.2)
            pc12(i,j) = 0;
        end
        if (pc13(i,j)>-0.2)&&(pc13(i,j)<0.2)
            pc13(i,j) = 0;
        end        
        if (pc14(i,j)>-0.2)&&(pc14(i,j)<0.2)
            pc14(i,j) = 0;
        end
        if (pc15(i,j)>-0.2)&&(pc15(i,j)<0.2)
            pc15(i,j) = 0;
        end
        if (pc16(i,j)>-0.2)&&(pc16(i,j)<0.2)
            pc16(i,j) = 0;
        end
    end
end
clear i j








