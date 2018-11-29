function te_result=transen(data_set,num_bin,max_tau)

%This function calculate the transfer entropy
%max_tau is maximum timelag.
%edited by Jewoo Hong(2011-04-11)
%edited by Juyeol Yun(2011-04-28)
[~,num_var]=size(data_set);

te_result=zeros(num_var,num_var,max_tau);

for i=1:num_var
    for j=1:num_var
        % H(Y(t-delta_t))
        pr_data=hist(data_set(1:end-1,j),num_bin)/sum(hist(data_set(1:end-1,j),num_bin));
        % H(Y(t),Y(t-delta_t))
        jpr_data=hist3(horzcat(data_set(2:end,j),data_set(1:end-1,j)),[num_bin,num_bin])...
            /sum(nansum(hist3(horzcat(data_set(2:end,j),data_set(1:end-1,j)),[num_bin,num_bin])));
        for tau=1:max_tau
            % H(X(t-tau*delta_t),Y(t-delta_t))
                xjpr_data=hist3(horzcat(data_set(1:end-tau,i),data_set(tau:end-1,j)),[num_bin,num_bin])...
                    /sum(nansum(hist3(horzcat(data_set(1:end-tau,i),data_set(tau:end-1,j)),[num_bin,num_bin])));
                tpr_data=hist4(horzcat(data_set(1:end-tau,i),data_set(tau+1:end,j),data_set(tau:end-1,j)),[num_bin,num_bin,num_bin])...
                    /sum(sum(nansum(hist4(horzcat(data_set(1:end-tau,i),data_set(tau+1:end,j),data_set(tau:end-1,j)),[num_bin,num_bin,num_bin]))));
            % H(X(t-tau*delta_t),Y(t),Y(t-delta_t))
            %result
            te_result(i,j,tau)=(-1)*nansum(xjpr_data(:).*log10(xjpr_data(:)))+(-1)*nansum(jpr_data(:).*log10(jpr_data(:)))-...
                (-1)*nansum(pr_data(:).*log10(pr_data(:)))-(-1)*nansum(tpr_data(:).*log10(tpr_data(:)));
        end
    end
end







