%% parallel
display('parallel')
matlabpool open 8;


%% name list
display('name list')

name_n = zeros(31,1);
temp = 0;
for i = 1:3
    temp=temp+1;
    name_n(temp,1) = 101991+i;
end
for i = 1:6
    temp=temp+1;
    name_n(temp,1) = 121993+i;
end
for i = 1:7
    temp=temp+1;
    name_n(temp,1) = 141996+i;
end
for i = 1:8
    temp=temp+1;
    name_n(temp,1) = 151999+i;
end
for i = 1:6
    temp=temp+1;
    name_n(temp,1) = 162003+i;
end
name_n(31,1) = 182010;
clear i temp


%% open
display('open')

for i = 1:31
    temp = ['open f' int2str(name_n(i,1)) '.mat;'];
    eval(temp);
    name = ['f' int2str(name_n(i,1)) ' = ans.f' int2str(name_n(i,1)) ';'];
    eval(name);
    clear ans
end

open('data_land.mat');
data_land = ans.data_land;

clear name temp i ans


%% land mask
display('land mask')

size_i = 1800;
size_j = 1800;
for i = 1:31
    temp = ['f' int2str(name_n(i,1)) '_land = zeros(size_i,size_j);'];
    eval(temp);
end
clear i temp
for i = 1:size_i
    for j = 1:size_j
        if data_land(i,j)>0
            for k = 1:31
                temp = ['f' int2str(name_n(k,1)) '_land(i,j) = f' int2str(name_n(k,1)) '(i,j);'];
                eval(temp);
            end
        end
    end
end
clear i j k temp


%% linear regression
display('linear regression')

pair = [3 4; 7 10; 8 11; 9 12; 13 17; 14 18; 15 19; 16 20; 21 25; 22 26; 23 27; 24 28];
pair_n = 12;
reg_coeff = zeros(pair_n,2);
reg_coeff2 = zeros(pair_n,1);
land_n = 0;

for i = 1:size_i
    for j = 1:size_j
        if data_land(i,j)==1
            land_n = land_n + 1;
        end
    end
end
clear i j
for i = 1:pair_n
    temp = ['regression_' int2str(i) ' = zeros(land_n,2);'];
    eval(temp);
end
clear i temp

land_n = 0;
for i = 1:size_i
    for j = 1:size_j
        if data_land(i,j)==1
            land_n = land_n + 1;
            for k = 1:pair_n
                temp = ['regression_' int2str(k) '(land_n,1) = f' int2str(name_n(pair(k,1),1)) '_land(i,j);'];
                eval(temp);
                temp = ['regression_' int2str(k) '(land_n,2) = f' int2str(name_n(pair(k,2),1)) '_land(i,j);'];
                eval(temp);
            end
        end        
    end
end
clear i j k temp

for i = 1:pair_n
    temp = ['polyfit(regression_' int2str(i) '(:,1), regression_' int2str(i) '(:,2),1);'];
    eval(temp);
    reg_coeff(i,1) = ans(1,1);
    reg_coeff(i,2) = ans(1,2);
end
clear i temp ans

for i = 1:pair_n
    xy = 0;
    x2 = 0;
    for j = 1:land_n
        temp = ['regression_' int2str(i) '(j,1)>0'];
        if eval(temp)
            temp = ['xy = xy + regression_' int2str(i) '(j,1)*regression_' int2str(i) '(j,2);'];
            eval(temp);
            temp = ['x2 = x2 + regression_' int2str(i) '(j,1)^2;'];
            eval(temp);
        end
    end
    reg_coeff2(i,1) = xy/x2;
end
clear i j xy x2
clear i temp pair pair_n land_n 


%% conversion
display('conversion')

for i = 1:31
    temp = ['f' int2str(name_n(i,1)) '_land_conv = f' int2str(name_n(i,1)) '_land;'];
    eval(temp);
end
clear i temp

slope = reg_coeff2(1,1)*reg_coeff2(2,1)*reg_coeff2(5,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f101992_land_conv(j,k) = f101992_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f101993_land_conv(j,k) = f101993_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f101994_land_conv(j,k) = f101994_land(j,k) * slope;
    end
end

slope = reg_coeff2(2,1)*reg_coeff2(5,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f121994_land_conv(j,k) = f121994_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f121995_land_conv(j,k) = f121995_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f121996_land_conv(j,k) = f121996_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f121997_land_conv(j,k) = f121997_land(j,k) * slope;
    end
end

slope = reg_coeff2(3,1)*reg_coeff2(5,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f121998_land_conv(j,k) = f121998_land(j,k) * slope;
    end
end

slope = reg_coeff2(4,1)*reg_coeff2(5,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f121999_land_conv(j,k) = f121999_land(j,k) * slope;
    end
end

slope = reg_coeff2(5,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f141997_land_conv(j,k) = f141997_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f141998_land_conv(j,k) = f141998_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f141999_land_conv(j,k) = f141999_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f142000_land_conv(j,k) = f142000_land(j,k) * slope;
    end
end

slope = reg_coeff2(6,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f142001_land_conv(j,k) = f142001_land(j,k) * slope;
    end
end

slope = reg_coeff2(7,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f142002_land_conv(j,k) = f142002_land(j,k) * slope;
    end
end

slope = reg_coeff2(8,1)*reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f142003_land_conv(j,k) = f142003_land(j,k) * slope;
    end
end

slope = reg_coeff2(9,1);
for j = 1:size_i
    parfor k = 1:size_j
        f152000_land_conv(j,k) = f152000_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f152001_land_conv(j,k) = f152001_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f152002_land_conv(j,k) = f152002_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f152003_land_conv(j,k) = f152003_land(j,k) * slope;
    end
end
for j = 1:size_i
    parfor k = 1:size_j
        f152004_land_conv(j,k) = f152004_land(j,k) * slope;
    end
end

slope = reg_coeff2(10,1);
for j = 1:size_i
    parfor k = 1:size_j
        f152005_land_conv(j,k) = f152005_land(j,k) * slope;
    end
end

slope = reg_coeff2(11,1);
for j = 1:size_i
    parfor k = 1:size_j
        f152006_land_conv(j,k) = f152006_land(j,k) * slope;
    end
end

slope = reg_coeff2(12,1);
for j = 1:size_i
    parfor k = 1:size_j
        f152007_land_conv(j,k) = f152007_land(j,k) * slope;
    end
end
clear i j k slope


%% year averaging
display('year averaging')

final_1992 = f101992_land_conv;
final_1993 = f101993_land_conv;
final_1995 = f121995_land_conv;
final_1996 = f121996_land_conv;
final_2008 = f162008_land_conv;
final_2009 = f162009_land_conv;

final_1994 = (f101994_land_conv + f121994_land_conv);
final_1994 = final_1994/2;

final_1997 = (f141997_land_conv + f121997_land_conv);
final_1997 = final_1997/2;
final_1998 = (f141998_land_conv + f121998_land_conv);
final_1998 = final_1998/2;
final_1999 = (f141999_land_conv + f121999_land_conv);
final_1999 = final_1999/2;

final_2000 = (f142000_land_conv + f152000_land_conv);
final_2000 = final_2000/2;
final_2001 = (f142001_land_conv + f152001_land_conv);
final_2001 = final_2001/2;
final_2002 = (f142002_land_conv + f152002_land_conv);
final_2002 = final_2002/2;
final_2003 = (f142003_land_conv + f152003_land_conv);
final_2003 = final_2003/2;

final_2004 = (f162004_land_conv + f152004_land_conv);
final_2004 = final_2004/2;
final_2005 = (f162005_land_conv + f152005_land_conv);
final_2005 = final_2005/2;
final_2006 = (f162006_land_conv + f152006_land_conv);
final_2006 = final_2006/2;
final_2007 = (f162007_land_conv + f152007_land_conv);
final_2007 = final_2007/2;


%% percentile conversion
display('percentile conversion')

for i = 1:size_i
    parfor j = 1:size_j
        if final_1992(i,j)>=63
            final_1992(i,j)=100;
        else
            final_1992(i,j) = final_1992(i,j)/63 *100;
        end
        if final_1993(i,j)>=63
            final_1993(i,j)=100;
        else
            final_1993(i,j) = final_1993(i,j)/63 *100;
        end
        if final_1994(i,j)>=63
            final_1994(i,j)=100;
        else
            final_1994(i,j) = final_1994(i,j)/63 *100;
        end
        if final_1995(i,j)>=63
            final_1995(i,j)=100;
        else
            final_1995(i,j) = final_1995(i,j)/63 *100;
        end
        if final_1996(i,j)>=63
            final_1996(i,j)=100;
        else
            final_1996(i,j) = final_1996(i,j)/63 *100;
        end
        if final_1997(i,j)>=63
            final_1997(i,j)=100;
        else
            final_1997(i,j) = final_1997(i,j)/63 *100;
        end
        if final_1998(i,j)>=63
            final_1998(i,j)=100;
        else
            final_1998(i,j) = final_1998(i,j)/63 *100;
        end        
        if final_1999(i,j)>=63
            final_1999(i,j)=100;
        else
            final_1999(i,j) = final_1999(i,j)/63 *100;
        end
        if final_2000(i,j)>=63
            final_2000(i,j)=100;
        else
            final_2000(i,j) = final_2000(i,j)/63 *100;
        end        
        if final_2001(i,j)>=63
            final_2001(i,j)=100;
        else
            final_2001(i,j) = final_2001(i,j)/63 *100;
        end                
        if final_2002(i,j)>=63
            final_2002(i,j)=100;
        else
            final_2002(i,j) = final_2002(i,j)/63 *100;
        end        
        if final_2003(i,j)>=63
            final_2003(i,j)=100;
        else
            final_2003(i,j) = final_2003(i,j)/63 *100;
        end        
        if final_2004(i,j)>=63
            final_2004(i,j)=100;
        else
            final_2004(i,j) = final_2004(i,j)/63 *100;
        end        
        if final_2005(i,j)>=63
            final_2005(i,j)=100;
        else
            final_2005(i,j) = final_2005(i,j)/63 *100;
        end        
        if final_2006(i,j)>=63
            final_2006(i,j)=100;
        else
            final_2006(i,j) = final_2006(i,j)/63 *100;
        end        
        if final_2007(i,j)>=63
            final_2007(i,j)=100;
        else
            final_2007(i,j) = final_2007(i,j)/63 *100;
        end        
        if final_2008(i,j)>=63
            final_2008(i,j)=100;
        else
            final_2008(i,j) = final_2008(i,j)/63 *100;
        end        
        if final_2009(i,j)>=63
            final_2009(i,j)=100;
        else
            final_2009(i,j) = final_2009(i,j)/63 *100;
        end                
    end
end
clear i j


%% test threshold
display('test threshold')

temp = final_1992;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1992 = land_n;

temp = final_1993;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1993 = land_n;


temp = final_1994;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1994 = land_n;

temp = final_1995;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1995 = land_n;

temp = final_1996;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1996 = land_n;

temp = final_1997;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1997 = land_n;

temp = final_1998;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1998 = land_n;

temp = final_1999;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th1999 = land_n;

temp = final_2000;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2000 = land_n;

temp = final_2001;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2001 = land_n;

temp = final_2002;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2002 = land_n;

temp = final_2003;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2003 = land_n;

temp = final_2004;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2004 = land_n;

temp = final_2005;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2005 = land_n;

temp = final_2006;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2006 = land_n;

temp = final_2007;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2007 = land_n;

temp = final_2008;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2008 = land_n;

temp = final_2009;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_i
        for j = 1:size_j
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
th2009 = land_n;

clear i j k temp land_n


%% KOREA separate
display('KOREA separate')

size_ki = 1420-775+1;
size_kj = 1320-540+1;

kor_1992 = zeros(size_ki,size_kj);
kor_1993 = zeros(size_ki,size_kj);
kor_1994 = zeros(size_ki,size_kj);
kor_1995 = zeros(size_ki,size_kj);
kor_1996 = zeros(size_ki,size_kj);
kor_1997 = zeros(size_ki,size_kj);
kor_1998 = zeros(size_ki,size_kj);
kor_1999 = zeros(size_ki,size_kj);
kor_2000 = zeros(size_ki,size_kj);
kor_2001 = zeros(size_ki,size_kj);
kor_2002 = zeros(size_ki,size_kj);
kor_2003 = zeros(size_ki,size_kj);
kor_2004 = zeros(size_ki,size_kj);
kor_2005 = zeros(size_ki,size_kj);
kor_2006 = zeros(size_ki,size_kj);
kor_2007 = zeros(size_ki,size_kj);
kor_2008 = zeros(size_ki,size_kj);
kor_2009 = zeros(size_ki,size_kj);

for i = 1:size_ki
    for j = 1:size_kj
        if ((i+775)>=1230)&&((j+540)>1090)
        elseif ((i+750)>=1360)&&((j+540)>865)
        else
            kor_1992(i,j) = final_1992(i+775+51,j+540+181);
            kor_1993(i,j) = final_1993(i+775+51,j+540+181);
            kor_1994(i,j) = final_1994(i+775+51,j+540+181);
            kor_1995(i,j) = final_1995(i+775+51,j+540+181);
            kor_1996(i,j) = final_1996(i+775+51,j+540+181);
            kor_1997(i,j) = final_1997(i+775+51,j+540+181);
            kor_1998(i,j) = final_1998(i+775+51,j+540+181);
            kor_1999(i,j) = final_1999(i+775+51,j+540+181);
            kor_2000(i,j) = final_2000(i+775+51,j+540+181);
            kor_2001(i,j) = final_2001(i+775+51,j+540+181);
            kor_2002(i,j) = final_2002(i+775+51,j+540+181);
            kor_2003(i,j) = final_2003(i+775+51,j+540+181);
            kor_2004(i,j) = final_2004(i+775+51,j+540+181);
            kor_2005(i,j) = final_2005(i+775+51,j+540+181);
            kor_2006(i,j) = final_2006(i+775+51,j+540+181);
            kor_2007(i,j) = final_2007(i+775+51,j+540+181);
            kor_2008(i,j) = final_2008(i+775+51,j+540+181);
            kor_2009(i,j) = final_2009(i+775+51,j+540+181);
        end
    end
end
clear i j 


%% test KOREA threshold
display('test KOREA threshold')

temp = kor_1992;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1992 = land_n;

temp = kor_1993;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1993 = land_n;

temp = kor_1994;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1994 = land_n;

temp = kor_1995;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1995 = land_n;

temp = kor_1996;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1996 = land_n;

temp = kor_1997;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1997 = land_n;

temp = kor_1998;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1998 = land_n;

temp = kor_1999;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1999 = land_n;

temp = kor_1995;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1995 = land_n;

temp = kor_1996;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1996 = land_n;

temp = kor_1997;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1997 = land_n;

temp = kor_1998;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1998 = land_n;

temp = kor_1999;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor1999 = land_n;

temp = kor_2000;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2000 = land_n;

temp = kor_2001;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2001 = land_n;

temp = kor_2002;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2002 = land_n;

temp = kor_2003;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2003 = land_n;

temp = kor_2004;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2004 = land_n;

temp = kor_2005;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2005 = land_n;

temp = kor_2006;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2006 = land_n;

temp = kor_2007;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2007 = land_n;

temp = kor_2008;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2008 = land_n;

temp = kor_2009;
land_n = zeros(100,1);
parfor k = 1:100
    for i = 1:size_ki
        for j = 1:size_kj
           if temp(i,j) > k
               land_n(k,1) = land_n(k,1) + 1;
           end
        end
    end
end
thkor2009 = land_n;

clear i j k temp land_n


%% perimeter of KOREA

temp = final_2009;
th = zeros(100,1);
for k = 1:100
    for i = 1:size_i
        for j = 1:size_j
            if temp(i,j)<k
                temp(i,j) = 0;
            end
        end
    end
    a=bwperim(temp);
    th(k,1) = sum(sum(a));
end

clear i j k a temp


% %% th ; 99
% display('th ; 99')
% 
% th = 99;
% kor_1992_th = kor_1992;
% kor_1993_th = kor_1993;
% kor_1994_th = kor_1994;
% kor_1995_th = kor_1995;
% kor_1996_th = kor_1996;
% kor_1997_th = kor_1997;
% kor_1998_th = kor_1998;
% kor_1999_th = kor_1999;
% kor_2000_th = kor_2000;
% kor_2001_th = kor_2001;
% kor_2002_th = kor_2002;
% kor_2003_th = kor_2003;
% kor_2004_th = kor_2004;
% kor_2005_th = kor_2005;
% kor_2006_th = kor_2006;
% kor_2007_th = kor_2007;
% kor_2008_th = kor_2008;
% kor_2009_th = kor_2009;
% 
% for i = 1:size_ki
%     parfor j = 1:size_kj
%         if kor_1992_th(i,j)<th
%             kor_1992_th(i,j) = 0;
%         end
%         if kor_1993_th(i,j)<th
%             kor_1993_th(i,j) = 0;
%         end
%         if kor_1994_th(i,j)<th
%             kor_1994_th(i,j) = 0;
%         end
%         if kor_1995_th(i,j)<th
%             kor_1995_th(i,j) = 0;
%         end
%         if kor_1996_th(i,j)<th
%             kor_1996_th(i,j) = 0;
%         end        
%         if kor_1997_th(i,j)<th
%             kor_1997_th(i,j) = 0;
%         end
%         if kor_1998_th(i,j)<th
%             kor_1998_th(i,j) = 0;
%         end
%         if kor_1999_th(i,j)<th
%             kor_1999_th(i,j) = 0;
%         end
%         if kor_2000_th(i,j)<th
%             kor_2000_th(i,j) = 0;
%         end
%         if kor_2001_th(i,j)<th
%             kor_2001_th(i,j) = 0;
%         end
%         if kor_2002_th(i,j)<th
%             kor_2002_th(i,j) = 0;
%         end
%         if kor_2003_th(i,j)<th
%             kor_2003_th(i,j) = 0;
%         end
%         if kor_2004_th(i,j)<th
%             kor_2004_th(i,j) = 0;
%         end
%         if kor_2005_th(i,j)<th
%             kor_2005_th(i,j) = 0;
%         end
%         if kor_2006_th(i,j)<th
%             kor_2006_th(i,j) = 0;
%         end
%         if kor_2007_th(i,j)<th
%             kor_2007_th(i,j) = 0;
%         end
%         if kor_2008_th(i,j)<th
%             kor_2008_th(i,j) = 0;
%         end
%         if kor_2009_th(i,j)<th
%             kor_2009_th(i,j) = 0;
%         end        
%     end
% end
% 
% clear i j th
% 

%% korea(seoul) mask data
display('korea(seoul) mask data')

data_marginal = bwperim(data_land);
data_mask = zeros(1800,1800);
for i = 1:size_i
    for j = 1:size_j
        data_mask(i,j) = data_marginal(i,j);
    end
end
clear data_marginal

for i = 1:6006
    lat_dec = mod(data_land_korea(i,1),1);
    lat_int = data_land_korea(i,1) - lat_dec;
    lon_dec = mod(data_land_korea(i,2),1);
    lon_int = data_land_korea(i,2) - lon_dec;
    
    lat_pos = int32((lat_int-45)*(-120)-lat_dec*120+mod(lat_dec*120,1));
    lon_pos = int32((lon_int-120)*120+1+(lon_dec*120-mod(lon_dec*120,1)));
    
    data_mask(lat_pos,lon_pos) = 1;
end

clear i lat_dec lat_int lon_dec lon_int lat_pos lon_pos


%% SEOUL separate
display('SEOUL separate')

size_si = 200;
size_sj = 200;

seoul_1992 = zeros(size_si,size_sj);
seoul_1993 = zeros(size_si,size_sj);
seoul_1994 = zeros(size_si,size_sj);
seoul_1995 = zeros(size_si,size_sj);
seoul_1996 = zeros(size_si,size_sj);
seoul_1997 = zeros(size_si,size_sj);
seoul_1998 = zeros(size_si,size_sj);
seoul_1999 = zeros(size_si,size_sj);
seoul_2000 = zeros(size_si,size_sj);
seoul_2001 = zeros(size_si,size_sj);
seoul_2002 = zeros(size_si,size_sj);
seoul_2003 = zeros(size_si,size_sj);
seoul_2004 = zeros(size_si,size_sj);
seoul_2005 = zeros(size_si,size_sj);
seoul_2006 = zeros(size_si,size_sj);
seoul_2007 = zeros(size_si,size_sj);
seoul_2008 = zeros(size_si,size_sj);
seoul_2009 = zeros(size_si,size_sj);
data_mask_seoul = zeros(size_si,size_sj);

for i = 1:size_si
    for j = 1:size_sj
        seoul_1992(i,j) = final_1992(i+775+51,j+540+181);
        seoul_1993(i,j) = final_1993(i+775+51,j+540+181);
        seoul_1994(i,j) = final_1994(i+775+51,j+540+181);
        seoul_1995(i,j) = final_1995(i+775+51,j+540+181);
        seoul_1996(i,j) = final_1996(i+775+51,j+540+181);
        seoul_1997(i,j) = final_1997(i+775+51,j+540+181);
        seoul_1998(i,j) = final_1998(i+775+51,j+540+181);
        seoul_1999(i,j) = final_1999(i+775+51,j+540+181);
        seoul_2000(i,j) = final_2000(i+775+51,j+540+181);
        seoul_2001(i,j) = final_2001(i+775+51,j+540+181);
        seoul_2002(i,j) = final_2002(i+775+51,j+540+181);
        seoul_2003(i,j) = final_2003(i+775+51,j+540+181);
        seoul_2004(i,j) = final_2004(i+775+51,j+540+181);
        seoul_2005(i,j) = final_2005(i+775+51,j+540+181);
        seoul_2006(i,j) = final_2006(i+775+51,j+540+181);
        seoul_2007(i,j) = final_2007(i+775+51,j+540+181);
        seoul_2008(i,j) = final_2008(i+775+51,j+540+181);
        seoul_2009(i,j) = final_2009(i+775+51,j+540+181);
        data_mask_seoul(i,j) = data_mask(i,j);
    end
end
clear i j 


%% exit parallel
display('exit parallel')
matlabpool close

