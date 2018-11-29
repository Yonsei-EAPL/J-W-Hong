TSTmean=zeros(4,1);
for i=1:4
    TSTmean(i,1)=mean(TSTm(i,:));
end

R_Rmean=zeros(4,1);
for i=1:4
    R_Rmean(i,1)=mean(R_R(i,:));
end

T_plusmean=zeros(4,1);
for i=1:4
    T_plusmean(i,1)=mean(T_plus(i,:));
end

T_minusmean=zeros(4,1);
for i=1:4
    T_minusmean(i,1)=mean(T_minus(i,:));
end

T_netmean=zeros(4,1);
for i=1:4
    T_netmean(i,1)=mean(T_net(i,:));
end

result=zeros(1,24);
for i=1:4
    result(1,i)=Hm(i,1);
end
for i=1:4
    result(1,i+4)=TSTmean(i,1);
end
for i=1:4
    result(1,i+8)=R_Rmean(i,1);
end
for i=1:4
    result(1,i+12)=T_plusmean(i,1);
end
for i=1:4
    result(1,i+16)=T_minusmean(i,1);
end
for i=1:4
    result(1,i+20)=T_netmean(i,1);
end

clear i A Ap Ap_rand H_ij H_ij_rand H_ji H_ji_rand H_result R R_R R_norm R_norm_rand R_rand TST TSTm T_minus T_net T_plus TSTmean R_Rmean T_plusmean T_minusmean T_netmean Hm graph w_te w_te_i