function z=montetrans(ano_data_set, i, j, tau, num_bin, num_ano_data)

% define
te_sample=0;

% transfer entropy
pr_tag_X=zeros(num_ano_data-tau,1);
pr_tag_Y=zeros(num_ano_data-tau,1);
pr_tag_Ym=zeros(num_ano_data-tau,1);
if i~=j
    % X
    X=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        X(k,1)=ano_data_set(k,i);
    end
    X=randsample(X,num_ano_data-tau);
    % Ym
    Ym=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        Ym(k,1)=ano_data_set(k+tau-1,j);
    end
    Ym=randsample(Ym,num_ano_data-tau);
    % Y
    Y=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        Y(k,1)=ano_data_set(k+tau,j);
    end
    Y=randsample(Y,num_ano_data-tau);
else
    data=randsample(ano_data_set(:,i),num_ano_data);
    % X
    X=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        X(k,1)=data(k,1);
    end
    % Ym
    Ym=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        Ym(k,1)=data(k+tau-1,1);
    end
    % Y
    Y=zeros(num_ano_data-tau,1);
    for k=1:num_ano_data-tau
        Y(k,1)=data(k+tau,1);
    end
end
% Y information
max_var_Y=max(Y);
for z=1:num_ano_data-tau
    if Y(z,1)==-99999
        Y(z,1)=99999;
    end
end
min_var_Y=min(Y);
for z=1:num_ano_data-tau
    if Y(z,1)==99999
        Y(z,1)=-99999;
    end
end   
dif_var_Y=max_var_Y-min_var_Y;
dif_cl_Y=dif_var_Y/num_bin;
% Ym information
max_var_Ym=max(Ym);
for z=1:num_ano_data-tau
    if Ym(z,1)==-99999
        Ym(z,1)=99999;
    end
end
min_var_Ym=min(Ym);
for z=1:num_ano_data-tau
    if Ym(z,1)==99999
        Ym(z,1)=-99999;
    end
end   
dif_var_Ym=max_var_Ym-min_var_Ym;
dif_cl_Ym=dif_var_Ym/num_bin;
% X information
max_var_X=max(X);
for z=1:num_ano_data-tau
    if X(z,1)==-99999
        X(z,1)=99999;
    end
end
min_var_X=min(X);
for z=1:num_ano_data-tau
    if X(z,1)==99999
        X(z,1)=-99999;
    end
end
dif_var_X=max_var_X-min_var_X;
dif_cl_X=dif_var_X/num_bin;
% p(Ym)
for k=1:num_ano_data-tau
    if Ym(k,1)~=-99999
        n=floor((Ym(k,1)-min_var_Ym)/dif_cl_Ym)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag_Ym(k,1)=n;
    end
end
pr_table_Ym=hist(pr_tag_Ym,0:1:num_bin);
num_Ym=sum(pr_table_Ym)-sum(pr_table_Ym(1,1));
pr_table_Ym=pr_table_Ym/num_Ym;
% p(Y,Ym)
for k=1:num_ano_data-tau
    if Y(k,1)~=-99999
        n=floor((Y(k,1)-min_var_Y)/dif_cl_Y)+1;
        if n==num_bin+1
            n=n-1;
        end
        % tag를 남겨서, 나중에 결합확률분포 계산 때, 다시 사용
        pr_tag_Y(k,1)=n;
    end
end
pr_table_Ycomb=hist3([pr_tag_Y pr_tag_Ym],{0:1:num_bin 0:1:num_bin});
num_Ycomb=sum(sum(pr_table_Ycomb(2:(num_bin+1),2:(num_bin+1))));
pr_table_Ycomb=pr_table_Ycomb/num_Ycomb;
%p(X,Y) and P(X,Y,Ym)
for k=1:num_ano_data-tau
    if X(k,1)~=-99999
        n=floor((X(k,1)-min_var_X)/dif_cl_X)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag_X(k,1)=n;
    end
end
pr_table_XYm=hist3([pr_tag_X pr_tag_Ym],{0:1:11 0:1:11});
num_XY=sum(sum(pr_table_XYm(2:(num_bin+1),2:(num_bin+1))));
pr_table_XYm=pr_table_XYm/num_XY;
pr_table_tri=histcnd([pr_tag_X pr_tag_Y pr_tag_Ym],{0:1:num_bin 0:1:num_bin 0:1:num_bin});
num_tri=sum(sum(sum(pr_table_tri(2:(num_bin+1),2:(num_bin+1),2:(num_bin+1)))));
pr_table_tri=pr_table_tri/num_tri;
% transfer entropy
for k=1:num_bin
    if pr_table_Ym(1,k+1)~=0
        te_sample=te_sample+pr_table_Ym(1,k+1)*log10(pr_table_Ym(1,k+1));
    end
    for s=1:num_bin
        if pr_table_Ycomb(k+1,s+1)~=0
            te_sample=te_sample-pr_table_Ycomb(k+1,s+1)*log10(pr_table_Ycomb(k+1,s+1));
        end
        if pr_table_XYm(s+1,k+1)~=0
            te_sample=te_sample-pr_table_XYm(s+1,k+1)*log10(pr_table_XYm(s+1,k+1));
        end
        for p=1:num_bin
            if pr_table_tri(s+1,p+1,k+1)~=0
                te_sample=te_sample+pr_table_tri(s+1,p+1,k+1)*log10(pr_table_tri(s+1,p+1,k+1));
            end
        end
    end
end

z=te_sample;