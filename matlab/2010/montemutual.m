function z=montemutual(ano_data_set, i, j, num_bin, num_ano_data)
% mutual information, monte carlo simulation 

% define
if i~=j
    data1=randsample(ano_data_set(:,i),num_ano_data);
    data2=randsample(ano_data_set(:,j),num_ano_data);
else
    data1=randsample(ano_data_set(:,i),num_ano_data);
    data2=data1;
end
pr_tag_X=zeros(num_ano_data,1);
pr_tag_Y=zeros(num_ano_data,1);
mi_sample=0;

% mutual information
% p(X)
max_var_X=max(data1);
for f=1:num_ano_data
    if data1(f,1)==-99999
        data1(f,1)=99999;
    end
end
min_var_X=min(data1);
for f=1:num_ano_data
    if data1(f,1)==99999
        data1(f,1)=-99999;
    end
end
dif_var_X=max_var_X-min_var_X;
dif_cl_X=dif_var_X/num_bin;
for m=1:num_ano_data
    if data1(m,1)~=-99999&&data2(m,1)~=-99999
        n=floor((data1(m,1)-min_var_X)/dif_cl_X)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag_X(m,1)=n;
    end
end
% p(Y), p(X,Y)
max_var_Y=max(data2);
for f=1:num_ano_data
    if data2(f,1)==-99999
        data2(f,1)=99999;
    end
end
min_var_Y=min(data2);
for f=1:num_ano_data
    if data2(f,1)==99999
        data2(f,1)=-99999;
    end
end
dif_var_Y=max_var_Y-min_var_Y;
dif_cl_Y=dif_var_Y/num_bin;
for m=1:num_ano_data
    if data2(m,1)~=-99999
        n=floor((data2(m,1)-min_var_Y)/dif_cl_Y)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag_Y(m,1)=n;
    end
end
pr_table_X=hist(pr_tag_X,0:1:num_bin);
num_X=sum(pr_table_X)-pr_table_X(1,1);
pr_table_X=pr_table_X/num_X;
pr_table_Y=hist(pr_tag_Y,0:1:num_bin);
num_Y=sum(pr_table_Y)-pr_table_Y(1,1);
pr_table_Y=pr_table_Y/num_Y;
pr_table_combine=hist3([pr_tag_X pr_tag_Y],{0:1:num_bin 0:1:num_bin});
num_comb=sum(sum(pr_table_combine))-sum(pr_table_combine(:,1))-sum(pr_table_combine(1,:))+pr_table_combine(1,1);
pr_table_combine=pr_table_combine/num_comb;

% mutual information
for f=1:num_bin
    for m=1:num_bin
        if pr_table_X(1,f+1)~=0&&pr_table_Y(1,m+1)~=0&&pr_table_combine(f+1,m+1)~=0
            mi_sample=mi_sample+pr_table_combine(f+1,m+1)*log10(pr_table_combine(f+1,m+1)/(pr_table_X(1,f+1)*pr_table_Y(1,m+1)));
        end               
    end
end

% resturn
z=mi_sample;