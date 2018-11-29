num_var=13;
clear AI AIr ATz Tau ano_data_set data_set mi_rand te_rand 

% H_result : normalized H
H_result=zeros(num_var,1);
for k=1:num_var
    H_result(k,1)=mi_result(k,k)/log10(num_bin);
end
clear k mi_result

% subgroup - 1 : turbulent, 2 : ABL, 3 : synoptic, 4 : total
% 1 : NEE, GPP, H, LE
% 2 : Rg, Precip
% 3 : RE, VPD, Tc, Ts, T, SWC, Pa
num_sys=3;

% using variable tag (tag)
tag=[3;1;1;3;1;2;2;3;3;3;3;1;3];

% Hm : mean H, [1,1]
% TST : sum A (A: adjacency metrix of TE), [tau,1]
% TSTm : mean A [tau,1]
% T_plus : mean gross information production (T[+]), [tau,1]
% T_minus : mean gross information consumtion (T[-]), [tau,1]
% T_net : mean net information production (T_net), [tau,1]
% Ap : A/TST, [num_var,num_var, tau]
% H_ij : conditional entropy, [tau,1]
% H_ji : conditional entropy, [tau,1]
% R : H_ij+H_ji, redundancy, [tau,1]
% R_norm : normalized R, [tau,1]
% R_norm_rand : 50times random sample

Hm=zeros(num_sys+1,1);
for i=1:num_sys+1
    if i==num_sys+1
        Hm(i,1)=sum(H_result)/num_var;
    else
        sum_H=0;
        sum_sys=0;
        for j=1:num_var
            if tag(j,1)==i
                sum_H=sum_H+H_result(j,1);
                sum_sys=sum_sys+1;
            end
        end
        Hm(i,1)=sum_H/sum_sys;
    end
end
clear i j sum_H sum_sys

TST=zeros(num_sys+1,max_tau);
TSTm=zeros(num_sys+1,max_tau);
% A : adjacency metrix of TE 
A=te_result/log10(num_bin);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau
            if A(i,j,k)<0
                A(i,j,k)=0;
            end
        end
    end
end
for i=1:num_sys+1
    for tau=1:max_tau
        l=0;
        m=0;
        if i==num_sys+1
            for j=1:num_var
                for k=1:num_var
                    l=l+A(j,k,tau);                   
                end
            end
            TST(i,tau)=l;
            TSTm(i,tau)=l/(num_var)^2;
        else
            for j=1:num_var
                if tag(j,1)==i
                    for k=1:num_var
                        if tag(k,1)==i
                            l=l+A(j,k,tau);
                            m=m+1;
                        end
                    end
                end
            end
            TST(i,tau)=l;
            TSTm(i,tau)=l/m^2;
        end
    end
end
clear i j k l m tau te_result

T_plus=zeros(num_sys+1,max_tau);
T_minus=zeros(num_sys+1,max_tau);
for i=1:num_sys+1
    for tau=1:max_tau
        l=0;
        m=0;
        if i==num_sys+1
            for j=1:num_var
                for k=1:num_var
                    l=l+A(j,k,tau);
                end
            end
            T_plus(i,tau)=l;
            T_minus(i,tau)=l;
        else
            for j=1:num_var
                if tag(j,1)==i
                    for k=1:num_var
                        l=l+A(j,k,tau);
                    end
                end
            end
            T_plus(i,tau)=l;
            for j=1:num_var
                for k=1:num_var
                    if tag(k,1)==i
                        m=m+A(j,k,tau);
                    end
                end
            end
            T_minus(i,tau)=m;
        end
    end
end
clear i j k l m n tau

T_net=T_plus-T_minus;


Ap=zeros(num_sys+1,num_var,num_var,max_tau);
for i=1:num_sys+1
    Ap(i,:,:,:)=A(1:num_var,1:num_var,:);
end
for i=1:num_sys+1
    if i==num_sys+1
    else
        for j=1:num_var
            if tag(j,1)==i
                for k=1:num_var
                    if tag(k,1)==i
                    else
                        Ap(i,j,k,:)=0;
                    end
                end
            else
                for k=1:num_var
                    Ap(i,j,k,:)=0;
                end
            end
        end
    end
end
for i=1:num_sys+1
    for j=1:max_tau
        Ap(i,:,:,j)=Ap(i,:,:,j)/TST(i,j);
    end
end
clear i j k

H_ij=zeros(num_sys+1,max_tau);
for i=1:num_sys+1
    for tau=1:max_tau
        l=0;
        if i==num_sys+1
            for j=1:num_var
                m=0;
                for k=1:num_var
                    m=m+Ap(i,j,k,tau);
                end
                for k=1:num_var
                    if m==0
                    else
                        l=l+Ap(i,j,k,tau)*log10(Ap(i,j,k,tau)/m);
                    end
                end
            end
            H_ij(i,tau)=-l;
        else
            for j=1:num_var
                m=0;
                if tag(j,1)==i
                    for k=1:num_var
                        if tag(k,1)==i
                            m=m+Ap(i,j,k,tau);
                        end
                    end
                    for k=1:num_var
                        if tag(k,1)==i
                            if m==0
                            else
                                l=l+Ap(i,j,k,tau)*log10(Ap(i,j,k,tau)/m);
                            end
                        end
                    end
                end
            end
            H_ij(i,tau)=-l;
        end
    end
end
clear i j k l m tau
H_ij=replace(H_ij,NaN,0);

H_ji=zeros(num_sys+1,max_tau); % j, k reverse from H_ij
for i=1:num_sys+1
    for tau=1:max_tau
        l=0;
        if i==num_sys+1
            for j=1:num_var
                m=0;
                for k=1:num_var
                    m=m+Ap(i,k,j,tau);
                end
                for k=1:num_var
                    if m==0
                    else
                        l=l+Ap(i,k,j,tau)*log10(Ap(i,k,j,tau)/m);
                    end
                end
            end
            H_ji(i,tau)=-l;
        else
            for j=1:num_var
                m=0;
                if tag(j,1)==i
                    for k=1:num_var
                        if tag(k,1)==i
                            m=m+Ap(i,k,j,tau);
                        end
                    end
                    for k=1:num_var
                        if tag(k,1)==i
                            if m==0
                            else
                                l=l+Ap(i,k,j,tau)*log10(Ap(i,k,j,tau)/m);
                            end
                        end
                    end
                end
            end
            H_ji(i,tau)=-l;
        end
    end
end
clear i j k l m tau
H_ji=replace(H_ji,NaN,0);

R=H_ij+H_ji;

R_norm=R;
for i=1:num_sys+1
    if i==num_sys+1
        R_norm(i,:)=R_norm(i,:)/(2*log10(num_var));
    else
        k=0;
        for j=1:num_var
            if tag(j,1)==i
                k=k+1;
            end
        end
        R_norm(i,:)=R_norm(i,:)/(2*log10(k));
    end
end
% clear R H_ij H_ji

% calculate R_rand(R_R) 50 times from random Ap 
R_R=zeros(num_sys+1,max_tau);
for o=1:50
    Ap_rand=zeros(num_sys+1,num_var,num_var,max_tau);
    for i=1:num_sys+1
        for j=1:max_tau
            k=Ap(i,:,:,j);
            l=k(randperm(num_var^2));
            for m=1:num_var
                for n=1:num_var
                    Ap_rand(i,m,n,j)=l(1,(m-1)*num_var+n);
                end
            end
        end
    end
    clear i j k l m n
    H_ij_rand=zeros(num_sys+1,max_tau);
    for i=1:num_sys+1
        for tau=1:max_tau
            l=0;
            if i==num_sys+1
                for j=1:num_var
                    m=0;
                    for k=1:num_var
                        m=m+Ap_rand(i,j,k,tau);
                    end
                    for k=1:num_var
                        if m==0
                        else
                            l=l+Ap_rand(i,j,k,tau)*log10(Ap_rand(i,j,k,tau)/m);
                        end
                    end
                end
                H_ij_rand(i,tau)=-l;
            else
                for j=1:num_var
                    m=0;
                    if tag(j,1)==i
                        for k=1:num_var
                            if tag(k,1)==i
                                m=m+Ap_rand(i,j,k,tau);
                            end
                        end
                        for k=1:num_var
                            if tag(k,1)==i
                                if m==0
                                else
                                    l=l+Ap_rand(i,j,k,tau)*log10(Ap_rand(i,j,k,tau)/m);
                                end
                            end
                        end
                    end
                end
                H_ij_rand(i,tau)=-l;
            end
        end
    end
    clear i j k l m tau
    H_ij_rand=replace(H_ij_rand,NaN,0);
    H_ji_rand=zeros(num_sys+1,max_tau); % j, k reverse from H_ij
    for i=1:num_sys+1
        for tau=1:max_tau
            l=0;
            if i==num_sys+1
                for j=1:num_var
                    m=0;
                    for k=1:num_var
                        m=m+Ap_rand(i,k,j,tau);
                    end
                    for k=1:num_var
                        if m==0
                        else
                            l=l+Ap_rand(i,k,j,tau)*log10(Ap_rand(i,k,j,tau)/m);
                        end
                    end
                end
                H_ji_rand(i,tau)=-l;
            else
                for j=1:num_var
                    m=0;
                    if tag(j,1)==i
                        for k=1:num_var
                            if tag(k,1)==i
                                m=m+Ap_rand(i,k,j,tau);
                            end
                        end
                        for k=1:num_var
                            if tag(k,1)==i
                                if m==0
                                else
                                    l=l+Ap_rand(i,k,j,tau)*log10(Ap_rand(i,k,j,tau)/m);
                                end
                            end
                        end
                    end
                end
                H_ji_rand(i,tau)=-l;
            end
        end
    end
    clear i j k l m tau
    H_ji_rand=replace(H_ji_rand,NaN,0);
    R_rand=H_ij_rand+H_ji_rand;
    R_norm_rand=R_rand;
    for i=1:num_sys+1
        if i==num_sys+1
            R_norm_rand(i,:)=R_norm_rand(i,:)/(2*log10(num_var));
        else
            k=0;
            for j=1:num_var
                if tag(j,1)==i
                    k=k+1;
                end
            end
            R_norm_rand(i,:)=R_norm_rand(i,:)/(2*log10(k));
        end
    end
    R_R=R_R+R_norm_rand;
end
clear i j k o num_sys
% clear R_norm_rand R_rand Ap_rand H_ij_rand H_ji_rand
R_R=R_R/50;
R_R=R_norm-R_R;
