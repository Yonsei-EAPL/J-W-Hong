% l = 21640;
% 
% info_ax = zeros(7,9);
% % title
% info_ax(1,1) = -2;
% info_ax(2,1) = 0;
% info_ax(3,1) = 0.1;
% info_ax(4,1) = 0.3;
% info_ax(5,1) = 0.5;
% info_ax(6,1) = 1;
% info_ax(7,1) = 2;
% % position_start&end
% for i = 1:l
%     temp= find(zol(:,1)==-2);
%     info_ax(1,2) = min(temp);
%     info_ax(1,3) = max(temp);
%     
%     temp= find(zol(:,1)==0);
%     info_ax(2,2) = min(temp);
%     info_ax(2,3) = max(temp);
%     
%     temp= find(zol(:,1)==0.1);
%     info_ax(3,2) = min(temp);
%     info_ax(3,3) = max(temp);
%     
%     temp= find(zol(:,1)==0.3);
%     info_ax(4,2) = min(temp);
%     info_ax(4,3) = max(temp);
%     
%     temp= find(zol(:,1)==0.5);
%     info_ax(5,2) = min(temp);
%     info_ax(5,3) = max(temp);
%     
%     temp= find(zol(:,1)==1);
%     info_ax(6,2) = min(temp);
%     info_ax(6,3) = max(temp);
%     
%     temp= find(zol(:,1)==2);
%     info_ax(7,2) = min(temp);
%     info_ax(7,3) = max(temp);
% end 
% clear i temp
% % ax_ start&end mean&min&max
% temp = ax(info_ax(1,2):info_ax(1,3),:);
% info_ax(1,4) = mean(temp(:,1));
% info_ax(1,5) = min(temp(:,1));
% info_ax(1,6) = max(temp(:,1));
% info_ax(1,7) = mean(temp(:,43));
% info_ax(1,8) = min(temp(:,43));
% info_ax(1,9) = max(temp(:,43));
% 
% temp = ax(info_ax(2,2):info_ax(2,3),:);
% info_ax(2,4) = mean(temp(:,1));
% info_ax(2,5) = min(temp(:,1));
% info_ax(2,6) = max(temp(:,1));
% info_ax(2,7) = mean(temp(:,43));
% info_ax(2,8) = min(temp(:,43));
% info_ax(2,9) = max(temp(:,43));
% 
% temp = ax(info_ax(3,2):info_ax(3,3),:);
% info_ax(3,4) = mean(temp(:,1));
% info_ax(3,5) = min(temp(:,1));
% info_ax(3,6) = max(temp(:,1));
% info_ax(3,7) = mean(temp(:,43));
% info_ax(3,8) = min(temp(:,43));
% info_ax(3,9) = max(temp(:,43));
% 
% temp = ax(info_ax(4,2):info_ax(4,3),:);
% info_ax(4,4) = mean(temp(:,1));
% info_ax(4,5) = min(temp(:,1));
% info_ax(4,6) = max(temp(:,1));
% info_ax(4,7) = mean(temp(:,43));
% info_ax(4,8) = min(temp(:,43));
% info_ax(4,9) = max(temp(:,43));
% 
% temp = ax(info_ax(5,2):info_ax(5,3),:);
% info_ax(5,4) = mean(temp(:,1));
% info_ax(5,5) = min(temp(:,1));
% info_ax(5,6) = max(temp(:,1));
% info_ax(5,7) = mean(temp(:,43));
% info_ax(5,8) = min(temp(:,43));
% info_ax(5,9) = max(temp(:,43));
% 
% temp = ax(info_ax(6,2):info_ax(6,3),:);
% info_ax(6,4) = mean(temp(:,1));
% info_ax(6,5) = min(temp(:,1));
% info_ax(6,6) = max(temp(:,1));
% info_ax(6,7) = mean(temp(:,43));
% info_ax(6,8) = min(temp(:,43));
% info_ax(6,9) = max(temp(:,43));
% 
% temp = ax(info_ax(7,2):info_ax(7,3),:);
% info_ax(7,4) = mean(temp(:,1));
% info_ax(7,5) = min(temp(:,1));
% info_ax(7,6) = max(temp(:,1));
% info_ax(7,7) = mean(temp(:,43));
% info_ax(7,8) = min(temp(:,43));
% info_ax(7,9) = max(temp(:,43));
% clear temp l

%% interpolation
ll = info_ax(1,3)-info_ax(1,2)+1;
fp_w_1 = zeros(ll, 41);
for i = 1:ll
    yi = interp1(ax(i,:),fp_w(i,:),ax_form1(1,:));
    for j = 1:41
        fp_w_1(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(2,3)-info_ax(2,2)+1;
fp_w_2 = zeros(ll, 27);
for i = 1:ll
    yi = interp1(ax(info_ax(2,2)-1+i,:),fp_w(info_ax(2,2)-1+i,:),ax_form2(1,:));
    for j = 1:27
        fp_w_2(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(3,3)-info_ax(3,2)+1;
fp_w_3 = zeros(ll, 49);
for i = 1:ll
    yi = interp1(ax(info_ax(3,2)-1+i,:),fp_w(info_ax(3,2)-1+i,:),ax_form3(1,:));
    for j = 1:49
        fp_w_3(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(4,3)-info_ax(4,2)+1;
fp_w_4 = zeros(ll, 52);
for i = 1:ll
    yi = interp1(ax(info_ax(4,2)-1+i,:),fp_w(info_ax(4,2)-1+i,:),ax_form4(1,:));
    for j = 1:52
        fp_w_4(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(5,3)-info_ax(5,2)+1;
fp_w_5 = zeros(ll, 54);
for i = 1:ll
    yi = interp1(ax(info_ax(5,2)-1+i,:),fp_w(info_ax(5,2)-1+i,:),ax_form5(1,:));
    for j = 1:54
        fp_w_5(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(6,3)-info_ax(6,2)+1;
fp_w_6 = zeros(ll, 50);
for i = 1:ll
    yi = interp1(ax(info_ax(6,2)-1+i,:),fp_w(info_ax(6,2)-1+i,:),ax_form6(1,:));
    for j = 1:50
        fp_w_6(i,j) = yi(1,j);
    end
end
clear i j yi ll
ll = info_ax(7,3)-info_ax(7,2)+1;
fp_w_7 = zeros(ll, 50);
for i = 1:ll
    yi = interp1(ax(info_ax(7,2)-1+i,:),fp_w(info_ax(7,2)-1+i,:),ax_form7(1,:));
    for j = 1:50
        fp_w_7(i,j) = yi(1,j);
    end
end
clear i j yi ll

%% arrange
fp_w2_mean = zeros(7,54);
fp_w2_err = zeros(7,54);

temp1 = mean(fp_w_1);
temp2 = std(fp_w_1);
l = length(fp_w_1);
for i = 1:41
    fp_w2_mean(1,i) = temp1(1,i);
    fp_w2_err(1,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_2);
temp2 = std(fp_w_2);
l = length(fp_w_2);
for i = 1:27
    fp_w2_mean(2,i) = temp1(1,i);
    fp_w2_err(2,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_3);
temp2 = std(fp_w_3);
l = length(fp_w_3);
for i = 1:49
    fp_w2_mean(3,i) = temp1(1,i);
    fp_w2_err(3,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_4);
temp2 = std(fp_w_4);
l = length(fp_w_4);
for i = 1:52
    fp_w2_mean(4,i) = temp1(1,i);
    fp_w2_err(4,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_5);
temp2 = std(fp_w_5);
l = length(fp_w_5);
for i = 1:54
    fp_w2_mean(5,i) = temp1(1,i);
    fp_w2_err(5,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_6);
temp2 = std(fp_w_6);
l = length(fp_w_6);
for i = 1:50
    fp_w2_mean(6,i) = temp1(1,i);
    fp_w2_err(6,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

temp1 = mean(fp_w_7);
temp2 = std(fp_w_7);
l = length(fp_w_7);
for i = 1:50
    fp_w2_mean(7,i) = temp1(1,i);
    fp_w2_err(7,i) = temp2(1,i)/sqrt(l);
end
clear i
clear temp1 temp2 l

clear fp_w_1 fp_w_2 fp_w_3 fp_w_4 fp_w_5 fp_w_6 fp_w_7


