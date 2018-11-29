% %%
% n_data = [1508568; 1687032; 1830840; 1940376; 2120016; 2190336; 2302896; 2444939; 2496648; 2543111; 2635183; 2608779; 2766000];
% n = sum(n_data);
% 
% data = zeros(n,7);
% j = 0;
% for i = 1:13
%     for k = 1:n_data(i,1)
%         j = j+1;
%         for l = 1:7
%             if i==1
%                 data(j,l) = x2001(k,l);
%             elseif i==2
%                 data(j,l) = x2002(k,l);
%             elseif i==3
%                 data(j,l) = x2003(k,l);
%             elseif i==4
%                 data(j,l) = x2004(k,l);
%             elseif i==5
%                 data(j,l) = x2005(k,l);
%             elseif i==6
%                 data(j,l) = x2006(k,l);
%             elseif i==7
%                 data(j,l) = x2007(k,l);
%             elseif i==8
%                 data(j,l) = x2008(k,l);
%             elseif i==9
%                 data(j,l) = x2009(k,l);
%             elseif i==10
%                 data(j,l) = x2010(k,l);
%             elseif i==11
%                 data(j,l) = x2011(k,l);
%             elseif i==12
%                 data(j,l) = x2012(k,l);
%             elseif i==13
%                 data(j,l) = x2013(k,l);
%             end    
%         end
%     end
% end
% clear i j k l

% %%
% n_stn = [111121;111122;111123;111124;111131;111141;111142;111151;111152;111153;111154;111161;111171;111181;111191;111201;111202;111212;111221;111231;111232;111241;111251;111261;111262;111263;111273;111274;111281;111291;111301;111311;131201;131202;131211];
% j= 0;
% n = 29074724;
% for i = 1:35
%     for k = 1:n
%         if data(k,1) == n_stn(i,1)
%             j = j+1;
%             n_stn(i,2) = n_stn(i,2)+1;
%         end
%     end
% end
% n = j;
% clear i j k
% 
% data2 = zeros(n,7);
% l = 0;
% for i = 1:35
%     for j = 1:29074724
%         if data(j,1) == n_stn(i,1)
%             l=l+1;
%             for k = 1:7
%                 data2(l,k) = data(j,k);
%             end
%         end
%     end
% end
% clear i j k l

%%
% time from excel (data2.mat)
i_stn = [111153,113952,1,1;111202,113952,1,1;111123,113952,2,1;111161,113952,2,1;111212,113952,2,1;111221,113952,2,1;111232,113952,2,1;111274,113952,2,1;111291,113952,2,1;111301,113952,2,1;111311,113952,2,1;111124,113952,3,1;111152,113952,3,1;111251,113952,3,1;111281,113952,3,1;111154,113952,4,1;111122,113952,5,1;111201,113952,6,1;111151,113952,7,2;111263,113952,7,2;111231,113952,8,2;111261,113952,8,2;111241,113952,9,2;111181,113952,9,2;111121,113952,10,2;131211,113952,11,2;111191,113952,12,3;111131,113952,13,3;111262,113952,14,3;131201,113952,15,3;111141,113952,16,3;111142,113952,16,4;131202,113952,17,4;111171,113952,18,4;111273,113952,19,4;];
n_data = 3988320;
n_lcz = 4;
n_stn = [18;8;5;4];

result = zeros(113952,20);
temp_sum = zeros(113952,20);
temp_n = zeros(113952,20);

for i = 1:35
    for j = 1:113952
        for k = 1:5
            if data2((i-1)*113952+j,k+2)>0
                temp_sum(j,(k-1)*4+i_stn(i,4)) = temp_sum(j,(k-1)*4+i_stn(i,4))+data2((i-1)*113952+j,k+2);
                temp_n(j,(k-1)*4+i_stn(i,4)) = temp_n(j,(k-1)*4+i_stn(i,4))+1;
            end
        end
    end
end
clear i j k
for i = 1:113952
    for j = 1:20
        if temp_n(i,j)>0
            result(i,j) = temp_sum(i,j)/temp_n(i,j);
        else
            result(i,j) = -999;
        end
    end
end
clear i j
for i = 1:113952
    for j = 1:4
        if result(i,j)>0
            result(i,j) = result(i,j)*1000;
        end
        if result(i,j+8)>0
            result(i,j+8) = result(i,j+8)*1000;
        end
        if result(i,j+12)>0
            result(i,j+12) = result(i,j+12)*1000;
        end        
    end
end
clear i j


% result = zeros(14,12,24,4,5);
% temp_sum1 = zeros(14,12,24,4);
% temp_n1 = zeros(14,12,24,4);
% temp_sum2 = zeros(14,12,24,4);
% temp_n2 = zeros(14,12,24,4);
% temp_sum3 = zeros(14,12,24,4);
% temp_n3 = zeros(14,12,24,4);
% temp_sum4 = zeros(14,12,24,4);
% temp_n4 = zeros(14,12,24,4);
% temp_sum5 = zeros(14,12,24,4);
% temp_n5 = zeros(14,12,24,4);

% 1: SO2 (ppm)
for i = 1:35
    for j = 1:113951 % time of 113952nd is 2014-01-01 00:00
        if data2((i-1)*113952+j,3)>0
            temp_sum1(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_sum1(time(j,1),time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,3);
            temp_n1(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_n1(time(j,1),time(j,2),time(j,3),i_stn(i,4))+1;
            temp_sum1(14,time(j,2),time(j,3),i_stn(i,4)) = temp_sum1(14,time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,3);
            temp_n1(14,time(j,2),time(j,3),i_stn(i,4)) = temp_n1(14,time(j,2),time(j,3),i_stn(i,4))+1;            
        end
    end
end
clear i j

% 2 : PM10 (ug m-3)
for i = 1:35
    for j = 1:113951 % time of 113952nd is 2014-01-01 00:00
        if data2((i-1)*113952+j,4)>0
            temp_sum2(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_sum2(time(j,1),time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,4);
            temp_n2(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_n2(time(j,1),time(j,2),time(j,3),i_stn(i,4))+1;
            temp_sum2(14,time(j,2),time(j,3),i_stn(i,4)) = temp_sum2(14,time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,4);
            temp_n2(14,time(j,2),time(j,3),i_stn(i,4)) = temp_n2(14,time(j,2),time(j,3),i_stn(i,4))+1;            
        end
    end
end
clear i j

% 3 : O3 (ppm)
for i = 1:35
    for j = 1:113951 % time of 113952nd is 2014-01-01 00:00
        if data2((i-1)*113952+j,5)>0
            temp_sum3(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_sum3(time(j,1),time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,5);
            temp_n3(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_n3(time(j,1),time(j,2),time(j,3),i_stn(i,4))+1;
            temp_sum3(14,time(j,2),time(j,3),i_stn(i,4)) = temp_sum3(14,time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,5);
            temp_n3(14,time(j,2),time(j,3),i_stn(i,4)) = temp_n3(14,time(j,2),time(j,3),i_stn(i,4))+1;            
        end
    end
end
clear i j

% 4 : NO2 (ppm)
for i = 1:35
    for j = 1:113951 % time of 113952nd is 2014-01-01 00:00
        if data2((i-1)*113952+j,6)>0
            temp_sum4(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_sum4(time(j,1),time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,6);
            temp_n4(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_n4(time(j,1),time(j,2),time(j,3),i_stn(i,4))+1;
            temp_sum4(14,time(j,2),time(j,3),i_stn(i,4)) = temp_sum4(14,time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,6);
            temp_n4(14,time(j,2),time(j,3),i_stn(i,4)) = temp_n4(14,time(j,2),time(j,3),i_stn(i,4))+1;                        
        end
    end
end
clear i j

% 5 : CO (ppm)
for i = 1:35
    for j = 1:113951 % time of 113952nd is 2014-01-01 00:00
        if data2((i-1)*113952+j,7)>0
            temp_sum5(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_sum5(time(j,1),time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,7);
            temp_n5(time(j,1),time(j,2),time(j,3),i_stn(i,4)) = temp_n5(time(j,1),time(j,2),time(j,3),i_stn(i,4))+1;
            temp_sum5(14,time(j,2),time(j,3),i_stn(i,4)) = temp_sum5(14,time(j,2),time(j,3),i_stn(i,4)) + data2((i-1)*113952+j,7);
            temp_n5(14,time(j,2),time(j,3),i_stn(i,4)) = temp_n5(14,time(j,2),time(j,3),i_stn(i,4))+1;            
        end
    end
end
clear i j

for i = 1:5
    for j = 1:14
        for k = 1:12
            for l = 1:24
                for m = 1:4
                    if i==1 %SO2
                        if temp_sum1(j,k,l,m)>0
                            result(j,k,l,m,i) = temp_sum1(j,k,l,m)/temp_n1(j,k,l,m);
                        else
                            result(j,k,l,m,i) = -999;
                        end
                    elseif i==2 %PM10
                        if temp_sum2(j,k,l,m)>0
                            result(j,k,l,m,i) = temp_sum2(j,k,l,m)/temp_n2(j,k,l,m);
                        else
                            result(j,k,l,m,i) = -999;
                        end                        
                    elseif i==3 %O3
                        if temp_sum3(j,k,l,m)>0
                            result(j,k,l,m,i) = temp_sum3(j,k,l,m)/temp_n3(j,k,l,m);
                        else
                            result(j,k,l,m,i) = -999;
                        end                        
                    elseif i==4 %NO2
                        if temp_sum4(j,k,l,m)>0
                            result(j,k,l,m,i) = temp_sum4(j,k,l,m)/temp_n4(j,k,l,m);
                        else
                            result(j,k,l,m,i) = -999;
                        end                        
                    elseif i==5 %CO
                        if temp_sum5(j,k,l,m)>0
                            result(j,k,l,m,i) = temp_sum5(j,k,l,m)/temp_n5(j,k,l,m);
                        else
                            result(j,k,l,m,i) = -999;
                        end                        
                    end
                end
            end
        end
    end
end
clear i j k l m



