function z=shannonentropy(ano_data_set, i, num_bin, num_ano_data);

%define
pr_table_X=zeros(num_bin,1);
H=0;

% p(X)
max_var_X=max(ano_data_set(:,i));
for f=1:num_ano_data
    if ano_data_set(f,i)==-99999
        ano_data_set(f,i)=99999;
    end
end
min_var_X=min(ano_data_set(:,i));
for f=1:num_ano_data
    if ano_data_set(f,i)==99999
        ano_data_set(f,i)=-99999;
    end
end
dif_var_X=max_var_X-min_var_X;
dif_cl_X=dif_var_X/num_bin;
for m=1:num_ano_data
    if ano_data_set(m,i)~=-99999
        n=floor((ano_data_set(m,i)-min_var_X)/dif_cl_X)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_table_X(n,1)=pr_table_X(n,1)+1;
    end
end
pr_table_X=pr_table_X/sum(pr_table_X);

%shannon entropy
for f=1:num_bin
    if pr_table_X(f,1)~=0
        H=H-pr_table_X(f,1)*log10(pr_table_X(f,1));
    end               
end

z=H;