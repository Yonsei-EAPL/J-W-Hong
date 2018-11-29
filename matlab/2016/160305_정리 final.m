% % % % l=14847;
% % % % for i = 1:l
% % % %     for j = 1:43
% % % %         x(i,j) = abs(x(i,j));
% % % %     end
% % % % end
% % % % clear i j
% % % 
% % % for i = 1:l
% % %     temp = sum(c(i,:));
% % %     for j = 1:43
% % %        c(i,j) =  c(i,j)/temp;
% % %     end
% % %     clear j
% % % 
% % %     temp = sum(t(i,:));
% % %     for j = 1:43
% % %        t(i,j) =  t(i,j)/temp;
% % %     end
% % %     clear j    
% % % 
% % %     temp = sum(u(i,:));
% % %     for j = 1:43
% % %        u(i,j) =  u(i,j)/temp;
% % %     end
% % %     clear j    
% % % 
% % %     temp = sum(v(i,:));
% % %     for j = 1:43
% % %        v(i,j) =  v(i,j)/temp;
% % %     end
% % %     clear j    
% % % 
% % %     temp = sum(w(i,:));
% % %     for j = 1:43
% % %        w(i,j) =  w(i,j)/temp;
% % %     end
% % %     clear j    
% % % 
% % %     temp = sum(x(i,:));
% % %     for j = 1:43
% % %        x(i,j) =  x(i,j)/temp;
% % %     end
% % %     clear j    
% % % end
% % % clear i temp
% % 
% % c2 = zeros(l,54);
% % t2 = zeros(l,54);
% % u2 = zeros(l,54);
% % v2 = zeros(l,54);
% % w2 = zeros(l,54);
% % x2 = zeros(l,54);
% % 
% % for i = 1:l
% %     temp = interp1(ax(i,:),c(i,:),ax2(1,:));
% %     for j = 1:54
% %         c2(i,j) = temp(1,j);
% %     end
% % 
% %     temp = interp1(ax(i,:),t(i,:),ax2(1,:));
% %     for j = 1:54
% %         t2(i,j) = temp(1,j);
% %     end    
% % 
% %     temp = interp1(ax(i,:),u(i,:),ax2(1,:));
% %     for j = 1:54
% %         u2(i,j) = temp(1,j);
% %     end    
% % 
% %     temp = interp1(ax(i,:),v(i,:),ax2(1,:));
% %     for j = 1:54
% %         v2(i,j) = temp(1,j);
% %     end   
% % 
% %     temp = interp1(ax(i,:),w(i,:),ax2(1,:));
% %     for j = 1:54
% %         w2(i,j) = temp(1,j);
% %     end    
% % 
% %     temp = interp1(ax(i,:),x(i,:),ax2(1,:));
% %     for j = 1:54
% %         x2(i,j) = temp(1,j);
% %     end    
% % end
% % clear i j temp
% 
% info_ax = zeros(7,3);
% % title
% info_ax(1,1) = -2;
% info_ax(2,1) = 0;
% info_ax(3,1) = 0.1;
% info_ax(4,1) = 0.3;
% info_ax(5,1) = 0.5;
% info_ax(6,1) = 1;
% info_ax(7,1) = 2;
% % position_start&end
% % for i = 1:l
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
% % end 
% clear i temp

% mean_c1 = mean(c2(info_ax(1,2):info_ax(1,3),:));
% mean_c2 = mean(c2(info_ax(2,2):info_ax(2,3),:));
% mean_c3 = mean(c2(info_ax(3,2):info_ax(3,3),:));
% mean_c4 = mean(c2(info_ax(4,2):info_ax(4,3),:));
% mean_c5 = mean(c2(info_ax(5,2):info_ax(5,3),:));
% mean_c6 = mean(c2(info_ax(6,2):info_ax(6,3),:));
% mean_c7 = mean(c2(info_ax(7,2):info_ax(7,3),:));
% mean_c = zeros(7,54);
% for i = 1:54
%     mean_c(1,i) = mean_c1(1,i);
%     mean_c(2,i) = mean_c2(1,i);
%     mean_c(3,i) = mean_c3(1,i);
%     mean_c(4,i) = mean_c4(1,i);
%     mean_c(5,i) = mean_c5(1,i);
%     mean_c(6,i) = mean_c6(1,i);
%     mean_c(7,i) = mean_c7(1,i);
% end
% clear i
% mean_c1 = std(c2(info_ax(1,2):info_ax(1,3),:));
% mean_c2 = std(c2(info_ax(2,2):info_ax(2,3),:));
% mean_c3 = std(c2(info_ax(3,2):info_ax(3,3),:));
% mean_c4 = std(c2(info_ax(4,2):info_ax(4,3),:));
% mean_c5 = std(c2(info_ax(5,2):info_ax(5,3),:));
% mean_c6 = std(c2(info_ax(6,2):info_ax(6,3),:));
% mean_c7 = std(c2(info_ax(7,2):info_ax(7,3),:));
% err_c = zeros(7,54);
% for i = 1:54
%     err_c(1,i) = mean_c1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
%     err_c(2,i) = mean_c2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
%     err_c(3,i) = mean_c3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
%     err_c(4,i) = mean_c4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
%     err_c(5,i) = mean_c5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
%     err_c(6,i) = mean_c6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
%     err_c(7,i) = mean_c7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
% end
% clear i mean_c1 mean_c2 mean_c3 mean_c4 mean_c5 mean_c6 mean_c7

mean_t1 = mean(t2(info_ax(1,2):info_ax(1,3),:));
mean_t2 = mean(t2(info_ax(2,2):info_ax(2,3),:));
mean_t3 = mean(t2(info_ax(3,2):info_ax(3,3),:));
mean_t4 = mean(t2(info_ax(4,2):info_ax(4,3),:));
mean_t5 = mean(t2(info_ax(5,2):info_ax(5,3),:));
mean_t6 = mean(t2(info_ax(6,2):info_ax(6,3),:));
mean_t7 = mean(t2(info_ax(7,2):info_ax(7,3),:));
mean_t = zeros(7,54);
for i = 1:54
    mean_t(1,i) = mean_t1(1,i);
    mean_t(2,i) = mean_t2(1,i);
    mean_t(3,i) = mean_t3(1,i);
    mean_t(4,i) = mean_t4(1,i);
    mean_t(5,i) = mean_t5(1,i);
    mean_t(6,i) = mean_t6(1,i);
    mean_t(7,i) = mean_t7(1,i);
end
clear i
mean_t1 = std(t2(info_ax(1,2):info_ax(1,3),:));
mean_t2 = std(t2(info_ax(2,2):info_ax(2,3),:));
mean_t3 = std(t2(info_ax(3,2):info_ax(3,3),:));
mean_t4 = std(t2(info_ax(4,2):info_ax(4,3),:));
mean_t5 = std(t2(info_ax(5,2):info_ax(5,3),:));
mean_t6 = std(t2(info_ax(6,2):info_ax(6,3),:));
mean_t7 = std(t2(info_ax(7,2):info_ax(7,3),:));
err_t = zeros(7,54);
for i = 1:54
    err_t(1,i) = mean_t1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
    err_t(2,i) = mean_t2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
    err_t(3,i) = mean_t3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
    err_t(4,i) = mean_t4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
    err_t(5,i) = mean_t5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
    err_t(6,i) = mean_t6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
    err_t(7,i) = mean_t7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
end
clear i mean_t1 mean_t2 mean_t3 mean_t4 mean_t5 mean_t6 mean_t7
mean_u1 = mean(u2(info_ax(1,2):info_ax(1,3),:));
mean_u2 = mean(u2(info_ax(2,2):info_ax(2,3),:));
mean_u3 = mean(u2(info_ax(3,2):info_ax(3,3),:));
mean_u4 = mean(u2(info_ax(4,2):info_ax(4,3),:));
mean_u5 = mean(u2(info_ax(5,2):info_ax(5,3),:));
mean_u6 = mean(u2(info_ax(6,2):info_ax(6,3),:));
mean_u7 = mean(u2(info_ax(7,2):info_ax(7,3),:));
mean_u = zeros(7,54);
for i = 1:54
    mean_u(1,i) = mean_u1(1,i);
    mean_u(2,i) = mean_u2(1,i);
    mean_u(3,i) = mean_u3(1,i);
    mean_u(4,i) = mean_u4(1,i);
    mean_u(5,i) = mean_u5(1,i);
    mean_u(6,i) = mean_u6(1,i);
    mean_u(7,i) = mean_u7(1,i);
end
clear i
mean_u1 = std(u2(info_ax(1,2):info_ax(1,3),:));
mean_u2 = std(u2(info_ax(2,2):info_ax(2,3),:));
mean_u3 = std(u2(info_ax(3,2):info_ax(3,3),:));
mean_u4 = std(u2(info_ax(4,2):info_ax(4,3),:));
mean_u5 = std(u2(info_ax(5,2):info_ax(5,3),:));
mean_u6 = std(u2(info_ax(6,2):info_ax(6,3),:));
mean_u7 = std(u2(info_ax(7,2):info_ax(7,3),:));
err_u = zeros(7,54);
for i = 1:54
    err_u(1,i) = mean_u1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
    err_u(2,i) = mean_u2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
    err_u(3,i) = mean_u3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
    err_u(4,i) = mean_u4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
    err_u(5,i) = mean_u5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
    err_u(6,i) = mean_u6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
    err_u(7,i) = mean_u7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
end
clear i mean_u1 mean_u2 mean_u3 mean_u4 mean_u5 mean_u6 mean_u7
mean_v1 = mean(v2(info_ax(1,2):info_ax(1,3),:));
mean_v2 = mean(v2(info_ax(2,2):info_ax(2,3),:));
mean_v3 = mean(v2(info_ax(3,2):info_ax(3,3),:));
mean_v4 = mean(v2(info_ax(4,2):info_ax(4,3),:));
mean_v5 = mean(v2(info_ax(5,2):info_ax(5,3),:));
mean_v6 = mean(v2(info_ax(6,2):info_ax(6,3),:));
mean_v7 = mean(v2(info_ax(7,2):info_ax(7,3),:));
mean_v = zeros(7,54);
for i = 1:54
    mean_v(1,i) = mean_v1(1,i);
    mean_v(2,i) = mean_v2(1,i);
    mean_v(3,i) = mean_v3(1,i);
    mean_v(4,i) = mean_v4(1,i);
    mean_v(5,i) = mean_v5(1,i);
    mean_v(6,i) = mean_v6(1,i);
    mean_v(7,i) = mean_v7(1,i);
end
clear i
mean_v1 = std(v2(info_ax(1,2):info_ax(1,3),:));
mean_v2 = std(v2(info_ax(2,2):info_ax(2,3),:));
mean_v3 = std(v2(info_ax(3,2):info_ax(3,3),:));
mean_v4 = std(v2(info_ax(4,2):info_ax(4,3),:));
mean_v5 = std(v2(info_ax(5,2):info_ax(5,3),:));
mean_v6 = std(v2(info_ax(6,2):info_ax(6,3),:));
mean_v7 = std(v2(info_ax(7,2):info_ax(7,3),:));
err_v = zeros(7,54);
for i = 1:54
    err_v(1,i) = mean_v1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
    err_v(2,i) = mean_v2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
    err_v(3,i) = mean_v3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
    err_v(4,i) = mean_v4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
    err_v(5,i) = mean_v5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
    err_v(6,i) = mean_v6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
    err_v(7,i) = mean_v7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
end
clear i mean_v1 mean_v2 mean_v3 mean_v4 mean_v5 mean_v6 mean_v7
mean_w1 = mean(w2(info_ax(1,2):info_ax(1,3),:));
mean_w2 = mean(w2(info_ax(2,2):info_ax(2,3),:));
mean_w3 = mean(w2(info_ax(3,2):info_ax(3,3),:));
mean_w4 = mean(w2(info_ax(4,2):info_ax(4,3),:));
mean_w5 = mean(w2(info_ax(5,2):info_ax(5,3),:));
mean_w6 = mean(w2(info_ax(6,2):info_ax(6,3),:));
mean_w7 = mean(w2(info_ax(7,2):info_ax(7,3),:));
mean_w = zeros(7,54);
for i = 1:54
    mean_w(1,i) = mean_w1(1,i);
    mean_w(2,i) = mean_w2(1,i);
    mean_w(3,i) = mean_w3(1,i);
    mean_w(4,i) = mean_w4(1,i);
    mean_w(5,i) = mean_w5(1,i);
    mean_w(6,i) = mean_w6(1,i);
    mean_w(7,i) = mean_w7(1,i);
end
clear i
mean_w1 = std(w2(info_ax(1,2):info_ax(1,3),:));
mean_w2 = std(w2(info_ax(2,2):info_ax(2,3),:));
mean_w3 = std(w2(info_ax(3,2):info_ax(3,3),:));
mean_w4 = std(w2(info_ax(4,2):info_ax(4,3),:));
mean_w5 = std(w2(info_ax(5,2):info_ax(5,3),:));
mean_w6 = std(w2(info_ax(6,2):info_ax(6,3),:));
mean_w7 = std(w2(info_ax(7,2):info_ax(7,3),:));
err_w = zeros(7,54);
for i = 1:54
    err_w(1,i) = mean_w1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
    err_w(2,i) = mean_w2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
    err_w(3,i) = mean_w3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
    err_w(4,i) = mean_w4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
    err_w(5,i) = mean_w5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
    err_w(6,i) = mean_w6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
    err_w(7,i) = mean_w7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
end
clear i mean_w1 mean_w2 mean_w3 mean_w4 mean_w5 mean_w6 mean_w7
mean_x1 = mean(x2(info_ax(1,2):info_ax(1,3),:));
mean_x2 = mean(x2(info_ax(2,2):info_ax(2,3),:));
mean_x3 = mean(x2(info_ax(3,2):info_ax(3,3),:));
mean_x4 = mean(x2(info_ax(4,2):info_ax(4,3),:));
mean_x5 = mean(x2(info_ax(5,2):info_ax(5,3),:));
mean_x6 = mean(x2(info_ax(6,2):info_ax(6,3),:));
mean_x7 = mean(x2(info_ax(7,2):info_ax(7,3),:));
mean_x = zeros(7,54);
for i = 1:54
    mean_x(1,i) = mean_x1(1,i);
    mean_x(2,i) = mean_x2(1,i);
    mean_x(3,i) = mean_x3(1,i);
    mean_x(4,i) = mean_x4(1,i);
    mean_x(5,i) = mean_x5(1,i);
    mean_x(6,i) = mean_x6(1,i);
    mean_x(7,i) = mean_x7(1,i);
end
clear i
mean_x1 = std(x2(info_ax(1,2):info_ax(1,3),:));
mean_x2 = std(x2(info_ax(2,2):info_ax(2,3),:));
mean_x3 = std(x2(info_ax(3,2):info_ax(3,3),:));
mean_x4 = std(x2(info_ax(4,2):info_ax(4,3),:));
mean_x5 = std(x2(info_ax(5,2):info_ax(5,3),:));
mean_x6 = std(x2(info_ax(6,2):info_ax(6,3),:));
mean_x7 = std(x2(info_ax(7,2):info_ax(7,3),:));
err_x = zeros(7,54);
for i = 1:54
    err_x(1,i) = mean_x1(1,i)/(info_ax(1,3)-info_ax(1,2)+1);
    err_x(2,i) = mean_x2(1,i)/(info_ax(2,3)-info_ax(2,2)+1);
    err_x(3,i) = mean_x3(1,i)/(info_ax(3,3)-info_ax(3,2)+1);
    err_x(4,i) = mean_x4(1,i)/(info_ax(4,3)-info_ax(4,2)+1);
    err_x(5,i) = mean_x5(1,i)/(info_ax(5,3)-info_ax(5,2)+1);
    err_x(6,i) = mean_x6(1,i)/(info_ax(6,3)-info_ax(6,2)+1);
    err_x(7,i) = mean_x7(1,i)/(info_ax(7,3)-info_ax(7,2)+1);
end
clear i mean_x1 mean_x2 mean_x3 mean_x4 mean_x5 mean_x6 mean_x7


