% y축 inverse해야 보는 방향과 남북방향이 일치함

temp=fxy_final;
sum(sum(temp))

temp2 = zeros(4001*4001,3);
for i = 1:4001
for j = 1:4001
temp2((i-1)*4001+j,1) = i;
temp2((i-1)*4001+j,2) = j;
temp2((i-1)*4001+j,3) = temp(i,j);
end
end
clear i j

temp2=sortrows(temp2,-3);
sum(temp2(:,3))

per = 0;
percen_loc=zeros(9,2);
for i = 1:16008001
per = per+temp2(i,3);
if percen_loc(1,1)==0
if per>10
percen_loc(1,1)=i;
percen_loc(1,2)=temp2(i,3);
end
end
if percen_loc(2,1)==0
if per>20
percen_loc(2,1)=i;
percen_loc(2,2)=temp2(i,3);
end
end
if percen_loc(3,1)==0
if per>30
percen_loc(3,1)=i;
percen_loc(3,2)=temp2(i,3);
end
end
if percen_loc(4,1)==0
if per>40
percen_loc(4,1)=i;
percen_loc(4,2)=temp2(i,3);
end
end
if percen_loc(5,1)==0
if per>50
percen_loc(5,1)=i;
percen_loc(5,2)=temp2(i,3);
end
end
if percen_loc(6,1)==0
if per>60
percen_loc(6,1)=i;
percen_loc(6,2)=temp2(i,3);
end
end
if percen_loc(7,1)==0
if per>70
percen_loc(7,1)=i;
percen_loc(7,2)=temp2(i,3);
end
end
if percen_loc(8,1)==0
if per>80
percen_loc(8,1)=i;
percen_loc(8,2)=temp2(i,3);
end
end
if percen_loc(9,1)==0
if per>90
percen_loc(9,1)=i;
percen_loc(9,2)=temp2(i,3);
end
end
end
clear i
clear i per

fp10=zeros(4001,4001);
fp20=zeros(4001,4001);
fp30=zeros(4001,4001);
fp40=zeros(4001,4001);
fp50=zeros(4001,4001);
fp60=zeros(4001,4001);
fp70=zeros(4001,4001);
fp80=zeros(4001,4001);
fp90=zeros(4001,4001);

per = 0;
for i = 1:16008001
per = per + temp2(i,3);
if per<10
fp10(temp2(i,1),temp2(i,2))=1;
fp20(temp2(i,1),temp2(i,2))=1;
fp30(temp2(i,1),temp2(i,2))=1;
fp40(temp2(i,1),temp2(i,2))=1;
fp50(temp2(i,1),temp2(i,2))=1;
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<20
fp20(temp2(i,1),temp2(i,2))=1;
fp30(temp2(i,1),temp2(i,2))=1;
fp40(temp2(i,1),temp2(i,2))=1;
fp50(temp2(i,1),temp2(i,2))=1;
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<30
fp30(temp2(i,1),temp2(i,2))=1;
fp40(temp2(i,1),temp2(i,2))=1;
fp50(temp2(i,1),temp2(i,2))=1;
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<40
fp40(temp2(i,1),temp2(i,2))=1;
fp50(temp2(i,1),temp2(i,2))=1;
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<50
fp50(temp2(i,1),temp2(i,2))=1;
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<60
fp60(temp2(i,1),temp2(i,2))=1;
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<70
fp70(temp2(i,1),temp2(i,2))=1;
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<80
fp80(temp2(i,1),temp2(i,2))=1;
fp90(temp2(i,1),temp2(i,2))=1;
elseif per<90
fp90(temp2(i,1),temp2(i,2))=1;
end
end
clear i



fp10_2 = bwperim(fp10);
fp20_2 = bwperim(fp20);
fp30_2 = bwperim(fp30);
fp40_2 = bwperim(fp40);
fp50_2 = bwperim(fp50);
fp60_2 = bwperim(fp60);
fp70_2 = bwperim(fp70);
fp80_2 = bwperim(fp80);
fp90_2 = bwperim(fp90);
figure()
spy(fp10_2);
hold on
spy(fp20_2);
spy(fp30_2);
spy(fp40_2);
spy(fp50_2);
spy(fp60_2);
spy(fp70_2);
spy(fp80_2);
spy(fp90_2);
