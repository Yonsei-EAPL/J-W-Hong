
%option 
num_data_day=48;
period_avg=5;
timestep=1800;
num_bin=11;
% data_set information
size_data=size(data_set);
num_var=size_data(1,2);
num_data=size_data(1,1);
clear size_data

% anomaly
ano_data_set=anomaly(data_set, num_var, num_data, period_avg, timestep);
clear period_avg timestep 
% anomaly information
size_data=size(ano_data_set);
num_ano_data=size_data(1,1);
num_day=num_ano_data(1,1)/num_data_day;
% �ڿ� ���ʿ��� ������ �����մϴ�.
clear size_data num_data_day

% mutual information
mi_result=zeros(num_var,num_var);
parfor i=1:num_var
    data1=num_ano_data(:,i);
    parfor j=1:num_var
        data2=num_ano_data(:,j);
        mi_result(i,j)=mutualinformation(data1, data2, num_bin, num_ano_data);
    end
end


% Part3.
% anomaly �����͸� �̿��Ͽ�, 
% ������ �ð��Ը�(max_tau) ������,
% �������� ���� �� �ִ� ��� ��쿡 ����, transfer entropy�� ����մϴ�.
% ����� te_result ��� �̸����� ��������� ���·� ��Ÿ���ϴ�.
% ���� �����ϰ��� �ϴ� �ð��Ը� �Ѱ��� "max_tau"�� �����մϴ�. 
% �⺻�� 36 (30�� ������ �� 18�ð�, ����&��)
max_tau=36;
% te ��꿡 ���Ǵ� �����
te_result=zeros(num_var,num_var,max_tau);
te_sample=0;
% Sink ���(Y)�� i�� �����Ͽ� ����
for i=1:num_var
    pr_t_comb_Y=zeros(num_bin,num_bin); % ���� ������ Y, Ym
    % Ym �� �ǹ̴� Y�� 1�� timestep ���� �����Ϳ��� ���մϴ�.(Y_t-dt)
    % ��꿡 �ʿ��� ǥ���� �ʱ�ȭ
    pr_t_Ym=zeros(num_bin,1);
    pr_t_comb_XY=zeros(num_bin,num_bin); %X, Y
    pr_t_tri=zeros(num_bin,num_bin,num_bin); %X, Y, Ym
    pr_tag_Y=zeros(num_ano_data-1,2); %Y, Ym
    % Y�� -1 timestep ���̰� ���� ������(Ym) ����
    % Ym�� ù��° �����ʹ� Ym�� t=2�� �� �����ͺ��� ���۵�.(t=1�� ����)
    Ym=zeros(num_ano_data-1,1);
    for j=1:num_ano_data-1
        Ym(j,1)=ano_data_set(j,i);
    end
    % Y ���� Ư�� �м�
    max_var_Y=max(ano_data_set(:,i));
    for z=1:num_ano_data
        if ano_data_set(z,i)==-99999
            ano_data_set(z,i)=99999;
        end
    end
    min_var_Y=min(ano_data_set(:,i));
    for z=1:num_ano_data
        if ano_data_set(z,i)==99999
            ano_data_set(z,i)=-99999;
        end
    end   
    dif_var_Y=max_var_Y-min_var_Y;
    dif_cl_Y=dif_var_Y/num_bin;
    % Ym Ư�� �м�
    max_var_Ym=max(Ym(:,1));
    for z=1:num_ano_data-1
        if Ym(z,1)==-99999
            Ym(z,1)=99999;
        end
    end
    min_var_Ym=min(Ym(:,1));
    for z=1:num_ano_data-1
        if Ym(z,1)==99999
            Ym(z,1)=-99999;
        end
    end   
    dif_var_Ym=max_var_Ym-min_var_Ym;
    dif_cl_Ym=dif_var_Ym/num_bin;
    % Ym �� Ȯ������ �ϼ�
    for j=1:num_ano_data-1
        if Ym(j,1)~=-99999
            n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
            if n==num_bin+1
                n=n-1;
            end
            pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
        % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
            pr_tag_Y(j,2)=n;
        end
    end
    pr_t_Ym=pr_t_Ym/sum(pr_t_Ym);
    % Y�� Ym�� ����Ȯ������ �ϼ�
    for j=2:num_ano_data
        if ano_data_set(j,i)~=-99999
            n=floor((ano_data_set(j,i)-min_var_Y)/dif_cl_Y)+1;
            if n==num_bin+1
                n=n-1;
            end
            % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
            pr_tag_Y(j-1,1)=n;
            if pr_tag_Y(j-1,2)~=0
                pr_t_comb_Y(n,pr_tag_Y(j-1,2))=pr_t_comb_Y(n,pr_tag_Y(j-1,2))+1;
            end
        end
    end
    % ������ ���� ������� 0���� ũ�� 1���� ���� Ȯ�� ���� �˴ϴ�.
    pr_t_comb_Y=pr_t_comb_Y/sum(sum(pr_t_comb_Y));
    % tau ��ŭ�� �ð����� ���� sink���(X) ~ m���� �� Ym�� ����Ȯ������ �ϼ�
    % source ���(X)�� m���� �����Ͽ� ����
    for m=1:num_var
        % �����ϴ� �ð��Ը� tau�� 1���� max_tau���� ����
        for tau=1:max_tau
            % X���� Ư�� ����
            Xt=zeros(num_ano_data-tau,1);
            for j=1:num_ano_data-tau
                Xt(j,1)=ano_data_set(j,m);
            end                
            max_var_X=max(Xt(:,1));
            for z=1:num_ano_data-tau
                if Xt(z,1)==-99999
                    Xt(z,1)=99999;
                end
            end
            min_var_X=min(Xt(:,1));
            for z=1:num_ano_data-tau
                if Xt(z,1)==99999
                    Xt(z,1)=-99999;
                end
            end
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            %XY����Ȯ�� �� XYYm����Ȯ�� �ϼ���Ű��
            for j=1:num_ano_data-tau
                if ano_data_set(j,m)~=-99999 
                    n=floor((ano_data_set(j,m)-min_var_X)/dif_cl_X)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    if pr_tag_Y(j+tau-1,2)~=0
                        pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))=pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))+1;
                        if pr_tag_Y(j+tau-1,1)~=0
                            pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))=pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))+1;                
                        end
                    end
                end
            end
            pr_t_comb_XY=pr_t_comb_XY/sum(sum(pr_t_comb_XY));
            pr_t_tri=pr_t_tri/sum(sum(sum(pr_t_tri)));
            % te ��� ����
            for j=1:num_bin
                if pr_t_Ym(j,1)~=0
                    te_sample=te_sample+pr_t_Ym(j,1)*log10(pr_t_Ym(j,1));
                end
                for s=1:num_bin
                    if pr_t_comb_XY(j,s)~=0
                        te_sample=te_sample-pr_t_comb_XY(j,s)*log10(pr_t_comb_XY(j,s));
                    end
                    if pr_t_comb_Y(j,s)~=0
                        te_sample=te_sample-pr_t_comb_Y(j,s)*log10(pr_t_comb_Y(j,s));
                    end
                    for p=1:num_bin
                        if pr_t_tri(j,s,p)~=0
                            te_sample=te_sample+pr_t_tri(j,s,p)*log10(pr_t_tri(j,s,p));
                        end
                    end
                end
            end
            % ����� �Է��ϰ�, ��꿡 �ݺ� ���Ǵ� ǥ �ʱ�ȭ
            te_result(m,i,tau)=te_sample;
            pr_t_tri=zeros(num_bin,num_bin,num_bin);
            pr_t_comb_XY=zeros(num_bin,num_bin);
            te_sample=0;     
        end
    end
end
% ���ʿ��� ������ �����մϴ�.
clear i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym Xt z
% ��.


% Part4. 
% Monte Carlo Simulator �κ� 
% randsample(); �Լ��� �̿��Ͽ�, 
% mutual information�� transfer entropy�� ���� �ݺ��Ͽ� �����ϰ�,
% ������ ���ڱ��� ����� �����Ͽ�,(�⺻���� 1000ȸ, cycle=1000;)
h=waitbar(0,'MI MonteCarlo, Please Wait...');
cycle=1000;
% �� ���� ���� mean �� standard deviation�� ����մϴ�.
% �����  �Ʒ��� ���� ������İ� ���� ���·� ��Ÿ���ϴ�.
mi_mean_rand=zeros(num_var,num_var); %��ü ���
mi_sigma_rand=zeros(num_var,num_var); %��ü ǥ������
mi_sample_rand=zeros(cycle,1); %������ mi ���
% Part2. �̿��ؼ�, mutual information, Monte Carlo simulation ����
% Ȯ������ ��� �� MI��꿡 �ʿ��� ���� ����
pr_table_X=zeros(num_bin,1);
pr_table_Y=zeros(num_bin,1);
pr_table_combine=zeros(num_bin,num_bin);
pr_tag=zeros(num_ano_data,1);
mi_sample=0;
X=zeros(num_ano_data,1);
Y=zeros(num_ano_data,1);
% Source ����(X)�� i�� �����Ͽ� ����
for i=1:num_var 
    X=ano_data_set(:,i);
    waitbar(i/num_var)
    for j=1:num_var
        Y=ano_data_set(:,j);
        for f=1:cycle %f�� ����Ƚ��
            rand_X=randsample(X,num_ano_data);        
            max_var_X=max(rand_X(:,1));
            for z=1:num_ano_data
                if rand_X(z,1)==-99999
                    rand_X(z,1)=99999;
                end
            end
            min_var_X=min(rand_X(:,1));
            for z=1:num_ano_data
                if rand_X(z,1)==99999
                    rand_X(z,1)=-99999;
                end
            end            
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            % ���� X�� Ȯ������ǥ�� �ϼ��ϰ�, ������ ������ ���� tag �ϼ�
            % tag�� �ڿ� ����Ȯ�������� ����� �� ���˴ϴ�.
            for m=1:num_ano_data
                if rand_X(m,1)~=-99999
                    n=floor((rand_X(m,1)-min_var_X)/dif_cl_X)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    pr_tag(m,1)=n;
                    pr_table_X(n,1)=pr_table_X(n,1)+1;
                end
            end
            pr_table_X=pr_table_X/sum(pr_table_X);
            % X�� Ȯ�������� �ϼ��Ǹ�, ������ sink ����(Y)������ MI�� ��� ����ϴ� ����
            rand_Y=randsample(Y,num_ano_data);
            max_var_Y=max(rand_Y(:,1));
            for z=1:num_ano_data
                if rand_Y(z,1)==-99999
                    rand_Y(z,1)=99999;
                end
            end
            min_var_Y=min(rand_Y(:,1));
            for z=1:num_ano_data
                if rand_Y(z,1)==99999
                    rand_Y(z,1)=-99999;
                end
            end    
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Y�� Ȯ�������� X,Y ����Ȯ������ �ϼ�
            for m=1:num_ano_data
                if rand_Y(m,1)~=-99999
                    n=floor((rand_Y(m,1)-min_var_Y)/dif_cl_Y)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    if pr_tag(m,1)~=0
                        pr_table_Y(n,1)=pr_table_Y(n,1)+1;
                        pr_table_combine(pr_tag(m,1),n)=pr_table_combine(pr_tag(m,1),n)+1;
                    end
                end
            end
            pr_table_Y=pr_table_Y/sum(pr_table_Y);
            pr_table_combine=pr_table_combine/sum(sum(pr_table_combine));
            % MI ��� ����
            for k=1:num_bin
                for m=1:num_bin
                    if pr_table_X(k,1)~=0&&pr_table_X(m,1)~=0&&pr_table_combine(k,m)~=0
                        mi_sample=mi_sample+pr_table_combine(k,m)*log10(pr_table_combine(k,m)/(pr_table_X(k,1)*pr_table_Y(m,1)));
                    end               
                end
            end
            % MI ������� ��� ǥ�� �Է��ϰ�, �ٽ� ����� ǥ���� �ʱ�ȭ
            mi_result_rand(f,1)=mi_sample;
            pr_table_Y=zeros(num_bin,1);
            pr_table_X=zeros(num_bin,1);
            pr_table_combine=zeros(num_bin,num_bin);
            pr_tag=zeros(num_ano_data,1);
            mi_sample=0;
            min_var_Y=0;
            max_var_Y=0;
            dif_var_Y=0;
            dif_cl_Y=0;
            min_var_X=0;
            max_var_X=0;
            dif_var_X=0;
            dif_cl_X=0;               
        end
        mi_mean_rand(i,j)=mean(mi_result_rand(:,1));
        mi_sigma_rand(i,j)=std(mi_result_rand(:,1));
        mi_sample_rand=zeros(cycle,1);
    end
end
close(h)
% ���ʿ��� ������ �����մϴ�.
clear pr_table_X pr_tag pr_table_Y pr_table_combine min_var_X min_var_Y max_var_X max_var_Y dif_var_X dif_var_Y dif_cl_X dif_cl_Y mi_sample i j k l m n X Y rand_X rand_Y mi_sample_rand mi_result_rand z
% Part2. �̿��Ͽ� mutual information�� Monte Carlo Simulation ��.
% Part3. �̿��Ͽ� transfer entropy�� Monte Carlo Simulation ����.
% ����� �Ʒ��� ���� ��������� ���·� ��Ÿ���ϴ�.
h=waitbar(0,'TE Monte Carlo, Please Wait...');
te_mean_rand=zeros(num_var,num_var,max_tau);
te_sigma_rand=zeros(num_var,num_var,max_tau);
te_sample_rand=zeros(cycle,max_tau);
% te ��꿡 ���Ǵ� �����
te_sample=0;
% Sink ���(Y)�� i�� �����Ͽ� ����
for i=1:num_var
    % source ���(X)�� m���� �����Ͽ� ����
    waitbar(i/num_var)
    for m=1:num_var
        X=ano_data_set(:,m);
        Y=ano_data_set(:,i);
        for f=1:cycle
            rand_X=randsample(X,num_ano_data);
            rand_Y=randsample(Y,num_ano_data);
            % Y�� 1 timestep ���̰� ���� ������(Ym) ����
            % Ym�� �ǹ̴� Y�� 1�� timestep ���� �����Ϳ��� ���մϴ�.(Y_t-dt) Ym�� t=2���� ����
            Ym=zeros(num_ano_data-1,1);
            for j=1:num_ano_data-1
                Ym(j,1)=rand_Y(j,1);
            end
            % ��꿡 �ʿ��� ǥ���� �ʱ�ȭ
            pr_t_comb_Y=zeros(num_bin,num_bin); % ���� ������ Y, Ym
            pr_t_Ym=zeros(num_bin,1);
            pr_t_comb_XY=zeros(num_bin,num_bin); %X, Y
            pr_t_tri=zeros(num_bin,num_bin,num_bin); %X, Y, Ym
            pr_tag_Y=zeros(num_ano_data-1,2); %Y, Ym            
            % Y ���� Ư�� �м�
            max_var_Y=max(rand_Y(:,1));
            for z=1:num_ano_data
                if rand_Y(z,1)==-99999
                    rand_Y(z,1)=99999;
                end
            end
            min_var_Y=min(rand_Y(:,1));
            for z=1:num_ano_data
                if rand_Y(z,1)==99999
                    rand_Y(z,1)=-99999;
                end
            end    
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Ym Ư�� �м�
            max_var_Ym=max(Ym(:,1));
            for z=1:num_ano_data-1
                if Ym(z,1)==-99999
                    Ym(z,1)=99999;
                end
            end
            min_var_Ym=min(Ym(:,1));
            for z=1:num_ano_data-1
                if Ym(z,1)==99999
                    Ym(z,1)=-99999;
                end
            end    
            dif_var_Ym=max_var_Y-min_var_Y;
            dif_cl_Ym=dif_var_Y/num_bin;
            % Ym �� Ȯ������ �ϼ�
            for j=1:num_ano_data-1
                if Ym(j,1)~=-99999
                    n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
                    % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
                    pr_tag_Y(j,2)=n;
                end
            end
            pr_t_Ym=pr_t_Ym/sum(pr_t_Ym);
            % Y�� Ym�� ����Ȯ������ �ϼ�
            for j=1:num_ano_data-1
                if rand_Y(j+1,1)~=-99999
                    n=floor((rand_Y(j+1,1)-min_var_Y)/dif_cl_Y)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
                    pr_tag_Y(j,1)=n;
                end
                if pr_tag_Y(j,2)~=0
                    pr_t_comb_Y(n,pr_tag_Y(j,2))=pr_t_comb_Y(n,pr_tag_Y(j,2))+1;
                end
            end
            % ������ ���� ������� 0���� ũ�� 1���� ���� Ȯ�� ���� �˴ϴ�.
            pr_t_comb_Y=pr_t_comb_Y/sum(sum(pr_t_comb_Y));
            % �����ϴ� �ð��Ը� tau�� 1���� max_tau���� ����
            for tau=1:max_tau
                % tau ��ŭ�� �ð����� ���� sink���(X) ~ m���� �� Ym�� ����Ȯ������ �ϼ�
                % X���� Ư�� ����
                Xt=zeros(num_ano_data-tau,1);
                for j=1:num_ano_data-tau
                    Xt(j,1)=ano_data_set(j,m);
                end                
                max_var_X=max(Xt(:,1));
                for z=1:num_ano_data-tau
                    if Xt(z,1)==-99999
                        Xt(z,1)=99999;
                    end
                end
                min_var_X=min(Xt(:,1));
                for z=1:num_ano_data-tau
                    if Xt(z,1)==99999
                        Xt(z,1)=-99999;
                    end
                end    
                dif_var_X=max_var_X-min_var_X;
                dif_cl_X=dif_var_X/num_bin;
                for j=1:num_ano_data-tau
                    if Xt(j,1)~=-99999
                        n=floor((Xt(j,1)-min_var_X)/dif_cl_X)+1;
                        if n==num_bin+1
                            n=n-1;
                        end
                        if pr_tag_Y(j+tau-1,2)~=0
                            pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))=pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))+1;
                            if pr_tag_Y(j+tau-1,1)~=0
                                pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))=pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))+1;                
                            end
                        end
                    end
                end
                pr_t_comb_XY=pr_t_comb_XY/sum(sum(pr_t_comb_XY));
                pr_t_tri=pr_t_tri/sum(sum(sum(pr_t_tri)));
                % te ��� ����
                for j=1:num_bin
                    if pr_t_Ym(j,1)~=0
                        te_sample=te_sample+pr_t_Ym(j,1)*log10(pr_t_Ym(j,1));
                    end
                    for s=1:num_bin
                        if pr_t_comb_XY(j,s)~=0
                            te_sample=te_sample-pr_t_comb_XY(j,s)*log10(pr_t_comb_XY(j,s));
                        end
                        if pr_t_comb_Y(j,s)~=0
                            te_sample=te_sample-pr_t_comb_Y(j,s)*log10(pr_t_comb_Y(j,s));
                        end
                        for p=1:num_bin
                            if pr_t_tri(j,s,p)~=0
                                te_sample=te_sample+pr_t_tri(j,s,p)*log10(pr_t_tri(j,s,p));
                            end
                        end
                    end
                end
                % ����� �Է��ϰ�, ��꿡 �ݺ� ���Ǵ� ǥ �ʱ�ȭ
                te_sample_rand(f,tau)=te_sample;
                pr_t_tri=zeros(num_bin,num_bin,num_bin);
                pr_t_comb_XY=zeros(num_bin,num_bin);
                te_sample=0;
            end
        end
        for j=1:max_tau
            te_mean_rand(m,i,j)=mean(te_sample_rand(:,j));
            te_sigma_rand(m,i,j)=std(te_sample_rand(:,j));
        end
        te_sample_rand=zeros(cycle,max_tau);
    end
end
% ���ʿ��� ������ �����մϴ�.
close(h)
clear f i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym X Y rand_X rand_Y te_sample_rand num_data_day num_day timestep period_avg Xt z h
% Part3. �̿� transfer entropy�� Monte Carlo Simulation ��.
% Part4. ��.


% Part5.
% Weighted Cut �κ� : ����� ���Ǽ�(mu+sigma)�� ���Ͽ�, ���� ��� ����
%   T > "mu(T_ss)+c*sigma(T_ss) = threshold of T"
%   95% : c=1.66
%   99% : c=2.36
% ���⼭�� �⺻�� c=1.66 (95%) ���
c=1.66;
% ����� mi_final(num_var,num_var) �� te_final(num_var,num_var,max_tau) �� ����
mi_final=mi_result;
te_final=te_result;
% ���Ǵ� �� ���� ����
mi_threshold=zeros(num_var,num_var);
te_threshold=zeros(num_var,num_var,max_tau);
% X(source)�� i, Y(sink)�� j�� rotation
for i=1:num_var
    for j=1:num_var
        mi_threshold(i,j)=mi_mean_rand(i,j)+c*mi_sigma_rand(i,j);
        for k=1:max_tau
            te_threshold(i,j,k)=te_mean_rand(i,j,k)+c*te_sigma_rand(i,j,k);    
        end
    end
end
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau
            if mi_final(i,j) < mi_threshold(i,j)
                mi_final(i,j)=0;
            end
            for k=1:max_tau
                if te_final(i,j,k) < te_threshold(i,j,k)
                    te_final(i,j,k)=0;
                end
            end
        end
    end
end
clear c i j k 
%clear mi_threshold te_threshold mi_mean_rand mi_sigma_rand te_mean_rand te_sigma_rand
% Part5. ��.


% Part6.
% normalizing �κ� : mi�� te�� max ���� log10(num_bin)���� ��ü�� ������
% ����� mi_final_nom�� te_final_nom ����
mi_final_nom=mi_final/log10(num_bin);
te_final_nom=te_final/log10(num_bin);
% Part6. ��.


% Part7.
% ������� ����� �κ�
% ����� AI(num_var,num_var)�� AIr(num_var,num_var) ATz(num_var,num_var), Tau(num_var,num_var), Tz(num_var,num_var,max_tau) ����
% ���� �������� ��Ʈ���� �ѷ��� ���մϴ�.
H_result=zeros(num_var,1);
pr_table_X=zeros(num_bin,1);
% Source ����(X)�� i�� �����Ͽ� ����
for i=1:num_var 
    H_X=0;
    max_var_X=max(ano_data_set(:,i));
    for z=1:num_ano_data
        if ano_data_set(z,i)==-99999
            ano_data_set(z,i)=99999;
        end
    end
    min_var_X=min(ano_data_set(:,i));
    for z=1:num_ano_data
        if ano_data_set(z,i)==99999
            ano_data_set(z,i)=-99999;
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
    for j=1:num_bin
        if pr_table_X(j,1)~=0
            H_X=H_X-pr_table_X(j,1)*log10(pr_table_X(j,1));
        end
    end
    H_result(i,1)=H_X;    
end
AI=zeros(i,j);
for i=1:num_var
    for j=1:num_var
           AI(i,j)=mi_result(i,j)/log10(num_bin)*100;
    end
end
AIr=zeros(i,j);
for i=1:num_var
    for j=1:num_var
        AIr(i,j)=100*mi_result(i,j)/H_result(j,1);
    end
end
Tz=zeros(num_var,num_var,max_tau);
Tz_final=zeros(num_var,num_var,max_tau);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau
            Tz(i,j,k)=te_result(i,j,k)/mi_result(i,j);
            if te_result(i,j,k)>te_threshold(i,j,k)
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
clear ans Tau_add i j k m n min_var_X max_var_X dif_var_X dif_cl_X pr_table_X H_X z
% Part7. ����.


% Part8. 
% �׷��� �ۼ��� �м��� ������, �����͸� �������ִ� �κ��Դϴ�.
% graph1 �� ����-������ ����� �ʱⰪ
% graph2 �� ����-�м��� threshold���� ����
% graph3 �� ����-������ weighted cut ������ ��..
% graph4 �� T/I�� ������ ��Ÿ����..
graph1=zeros(max_tau+1,num_var^2);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau+1
            if k==1
                graph1(1,(i-1)*num_var+j)=mi_result(i,j);
            else
                graph1(k,(i-1)*num_var+j)=te_result(i,j,k-1);
            end
        end
    end
end
graph2=zeros(max_tau+1,num_var^2);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau+1
            if k==1
                graph2(1,(i-1)*num_var+j)=mi_threshold(i,j);
            else
                graph2(k,(i-1)*num_var+j)=te_threshold(i,j,k-1);
            end
        end
    end
end
graph3=zeros(max_tau+1,num_var^2);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau+1
            if k==1
                graph3(1,(i-1)*num_var+j)=mi_final(i,j);
            else
                graph3(k,(i-1)*num_var+j)=te_final(i,j,k-1);
            end
        end
    end
end
graph4=zeros(max_tau+1,num_var^2);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau+1
            if k==1
                graph4(1,(i-1)*num_var+j)=mi_result(i,j)/mi_result(i,j);
            else
                graph4(k,(i-1)*num_var+j)=te_result(i,j,k-1)/mi_result(i,j);
            end
        end
    end
end
graph5=zeros(max_tau+1,num_var^2);
for i=1:num_var
    for j=1:num_var
        for k=1:max_tau+1
            if k==1
                graph5(1,(i-1)*num_var+j)=mi_threshold(i,j)/mi_result(i,j);
            else
                graph5(k,(i-1)*num_var+j)=te_threshold(i,j,k-1)/mi_result(i,j);
            end
        end
    end
end
clear i j k
%Part8. ��. 


% ���� ����� ���̻� �ʿ���ٸ� �Ʒ��κб��� �����ϸ�, ����� �����ϴ�.
clear num_data_day num_day timestep period_avg
%clear num_var num_ano_data num_bin max_tau