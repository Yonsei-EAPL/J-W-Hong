%after estimation
clear confidence iteration nom_w_mi nom_w_te nom_w_te_i period timestep w_rmi te_result  erf_te_rand
te_result=te(:,:,2:37);
clear te
mi_rand=erf_mi_rand;
clear erf_mi_rand
te_rand=erf_te(:,:,2:37);
clear erf_te

%adjacency matrix
AI=zeros(num_var,num_var);
for i=1:num_var
    for j=1:num_var
           AI(i,j)=mi_result(i,j)/log10(num_bin)*100;
    end
end
AIr=zeros(num_var,num_var);
for i=1:num_var
    for j=1:num_var
        AIr(i,j)=100*mi_result(i,j)/mi_result(j,j);
    end
end

Tz=zeros(num_var,num_var,max_tau);
Tz_final=zeros(num_var,num_var,max_tau);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau
            Tz(i,j,k)=te_result(i,j,k)/mi_result(i,j);
            if te_result(i,j,k)>te_rand(i,j,k)
                Tz_final(i,j,k)=Tz(i,j,k);
            else
                Tz_final(i,j,k)=0;
            end
        end
    end
end
ATz=zeros(num_var,num_var);
for i=1:num_var
    for j=1:num_var
        ATz(i,j)=max(Tz_final(i,j,:));
    end
end
Tau=zeros(i,j,4);
for i=1:num_var
    for j=1:num_var
        Tau_add=find(Tz_final(i,j,:)==max(Tz_final(i,j,:)));
        size(Tau_add);
        if ans(1,1)==36
            Tau(i,j,:)=0;
        else
            Tau(i,j,4)=Tau_add;        
            Tau_add=find(Tz_final(i,j,:)>0);
            Tau(i,j,1)=min(Tau_add);
            Tau(i,j,2)=max(Tau_add);
            size(Tau_add);
            Tau(i,j,3)=ans(1,1);
        end
    end
end
clear ans Tau_add i j Tz Tz_final k