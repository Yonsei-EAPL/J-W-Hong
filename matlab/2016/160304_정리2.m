%% interpolation
ll = info_ax(1,3)-info_ax(1,2)+1;
fp_coab_1 = zeros(ll, 41);
for i = 1:ll
    yi = interp1(ax(i,:),fp_coab(i,:),ax_form1(1,:));
    for j = 1:41
        fp_coab_1(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(2,3)-info_ax(2,2)+1;
fp_coab_2 = zeros(ll, 27);
for i = 1:ll
    yi = interp1(ax(info_ax(2,2)-1+i,:),fp_coab(info_ax(2,2)-1+i,:),ax_form2(1,:));
    for j = 1:27
        fp_coab_2(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(3,3)-info_ax(3,2)+1;
fp_coab_3 = zeros(ll, 49);
for i = 1:ll
    yi = interp1(ax(info_ax(3,2)-1+i,:),fp_coab(info_ax(3,2)-1+i,:),ax_form3(1,:));
    for j = 1:49
        fp_coab_3(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(4,3)-info_ax(4,2)+1;
fp_coab_4 = zeros(ll, 52);
for i = 1:ll
    yi = interp1(ax(info_ax(4,2)-1+i,:),fp_coab(info_ax(4,2)-1+i,:),ax_form4(1,:));
    for j = 1:52
        fp_coab_4(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(5,3)-info_ax(5,2)+1;
fp_coab_5 = zeros(ll, 54);
for i = 1:ll
    yi = interp1(ax(info_ax(5,2)-1+i,:),fp_coab(info_ax(5,2)-1+i,:),ax_form5(1,:));
    for j = 1:54
        fp_coab_5(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(6,3)-info_ax(6,2)+1;
fp_coab_6 = zeros(ll, 50);
for i = 1:ll
    yi = interp1(ax(info_ax(6,2)-1+i,:),fp_coab(info_ax(6,2)-1+i,:),ax_form6(1,:));
    for j = 1:50
        fp_coab_6(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(7,3)-info_ax(7,2)+1;
fp_coab_7 = zeros(ll, 50);
for i = 1:ll
    yi = interp1(ax(info_ax(7,2)-1+i,:),fp_coab(info_ax(7,2)-1+i,:),ax_form7(1,:));
    for j = 1:50
        fp_coab_7(i,j) = yi(1,j);
    end
end
clear i j yi ll

%% arrange
fp_coab2_mean = zeros(7,54);
fp_coab2_err = zeros(7,54);

temp1 = mean(fp_coab_1);
temp2 = std(fp_coab_1);
l = length(fp_coab_1);
for i = 1:41
    fp_coab2_mean(1,i) = temp1(1,i);
    fp_coab2_err(1,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_2);
temp2 = std(fp_coab_2);
l = length(fp_coab_2);
for i = 1:27
    fp_coab2_mean(2,i) = temp1(1,i);
    fp_coab2_err(2,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_3);
temp2 = std(fp_coab_3);
l = length(fp_coab_3);
for i = 1:49
    fp_coab2_mean(3,i) = temp1(1,i);
    fp_coab2_err(3,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_4);
temp2 = std(fp_coab_4);
l = length(fp_coab_4);
for i = 1:52
    fp_coab2_mean(4,i) = temp1(1,i);
    fp_coab2_err(4,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_5);
temp2 = std(fp_coab_5);
l = length(fp_coab_5);
for i = 1:54
    fp_coab2_mean(5,i) = temp1(1,i);
    fp_coab2_err(5,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_6);
temp2 = std(fp_coab_6);
l = length(fp_coab_6);
for i = 1:50
    fp_coab2_mean(6,i) = temp1(1,i);
    fp_coab2_err(6,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_coab_7);
temp2 = std(fp_coab_7);
l = length(fp_coab_7);
for i = 1:50
    fp_coab2_mean(7,i) = temp1(1,i);
    fp_coab2_err(7,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

clear fp_coab_1 fp_coab_2 fp_coab_3 fp_coab_4 fp_coab_5 fp_coab_6 fp_coab_7
