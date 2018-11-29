% �� ���α׷��� �ۼ��� ���������� �Դϴ�. (made by Je-woo Hong)
% �ǰ� �ֽø� �ݿ��Ͽ� ������ ����ڽ��ϴ�.

% Version 1.0 (2011-03-09, Je-woo Hong)
% ����� data�� matlab ���ο��� �����Ͽ� �Űܿͼ� ����ؾ� �մϴ�.
% ��� ���鵵, matlab ���ο� ��µ˴ϴ�.
% �������� ���׷��̵� ����
% - data �� result�� ����(�Ǵ� notepad)���� �а� ����
% - �ڵ����� �Ѵ޾� ��� ����ϱ�
% - ����&�� 2��° ���� ������ ���
% - ����&�� 1��° ���� ����ī���� �ùķ��̼� �� ������� �� �߰����� �߰�(Ver1.1���� �ذ�)

% Version 1.1 (2011-04-28, Je-woo Hong)
% �߰� �� ���� �� ����
% - Part4. MonteCarlo Simulation �߰�
% - Part3. Transfer Entropy ��� �� ����(data �ð�����) ����
% - Part5. Weighted Cut �߰�
% - Part6. Normalizing �߰�
% - Part7. ������� �߰�


% Part1.
% data_set�� anomaly�� ����� �κ��Դϴ�.
% ������� data_set�� �����ϰ�, data�� �Է��ϼ���.
% ex) data_set=zeros(�����ͼ�, ������);
% ���� anomaly�� ����� ���� ����ϴ� �Ⱓ�� �����մϴ�. "period_avg"
% �⺻���� 5�� (5��(���ð� ���� 4��)���� �����ð��� �����͸� ��ճ��� ���)
period_avg=5;
% �������� "timestep" ���� (���� : sec), �⺻���� 30min = 1800sec
timestep=1800;
% �������� �⺻ ���� ���� �����մϴ�.
size_data=size(data_set);
num_var=size_data(1,2)
num_data=size_data(1,1)
num_data_day=24*3600/timestep
num_day=num_data(1,1)/num_data_day
avg_sample=0;
ano_data_set=zeros(num_data-num_data_day*(period_avg-1),num_var);
% anomaly �����͸� �����մϴ�.
for i=1:num_var
    for j=(period_avg-1)*num_data_day+1:num_data
        for k=1:period_avg
            avg_sample=avg_sample+data_set(j-(k-1)*num_data_day,i);
        end
        avg_sample=avg_sample/period_avg;
        ano_data_set(j-(period_avg-1)*num_data_day,i)=data_set(j,i)-avg_sample;
        avg_sample=0;
    end
end
% anomaly ������ ������ �����մϴ�.
size_data=size(ano_data_set);
num_ano_data=size_data(1,1)
num_day=num_ano_data(1,1)/num_data_day
% �ڿ� ���ʿ��� ������ �����մϴ�.
clear i j k avg_sample size_data num_data
% Part1. ��.


% Part2.
% anomaly �����͸� �̿��Ͽ�, 
% �������� ���� �� �ִ� ��� �����, mutual information�� ����մϴ�.
% ����� mi_result ��� �̸����� ��������� ���·� ��Ÿ���ϴ�.
% ���� bin ���� "num_bin"�� �������ּ���. �⺻���� 11 (����&��)
num_bin=11;
% Ȯ������ ��� �� MI��꿡 �ʿ��� ���� ����
pr_table_X=zeros(num_bin,1);
pr_table_Y=zeros(num_bin,1);
pr_table_combine=zeros(num_bin,num_bin);
pr_tag=zeros(num_ano_data,1);
mi_result=zeros(num_var,num_var);
mi_sample=0;
% Source ����(X)�� i�� �����Ͽ� ����
for i=1:num_var 
    min_var_X=min(ano_data_set(:,i));
    max_var_X=max(ano_data_set(:,i));
    dif_var_X=max_var_X-min_var_X;
    dif_cl_X=dif_var_X/num_bin;
    % ���� X�� Ȯ������ǥ�� �ϼ��ϰ�, ������ ������ ���� tag �ϼ�
    % tag�� �ڿ� ����Ȯ�������� ����� �� ���˴ϴ�.
    for m=1:num_ano_data
        n=floor((ano_data_set(m,i)-min_var_X)/dif_cl_X)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag(m,1)=n;
        pr_table_X(n,1)=pr_table_X(n,1)+1;
    end
    pr_table_X=pr_table_X/num_ano_data;
    % X�� Ȯ�������� �ϼ��Ǹ�, ������ sink ����(Y)������ MI�� ��� ����ϴ� ����
    for j=1:num_var
        min_var_Y=min(ano_data_set(:,j));
        max_var_Y=max(ano_data_set(:,j));
        dif_var_Y=max_var_Y-min_var_Y;
        dif_cl_Y=dif_var_Y/num_bin;
        % Y�� Ȯ�������� X,Y ����Ȯ������ �ϼ�
        for l=1:num_ano_data
            n=floor((ano_data_set(l,j)-min_var_Y)/dif_cl_Y)+1;
            if n==num_bin+1
                n=n-1;
            end
            pr_table_Y(n,1)=pr_table_Y(n,1)+1;
            pr_table_combine(pr_tag(l,1),n)=pr_table_combine(pr_tag(l,1),n)+1;
        end
        pr_table_Y=pr_table_Y/num_ano_data;
        pr_table_combine=pr_table_combine/num_ano_data;
        % MI ��� ����
        for k=1:num_bin
            for m=1:num_bin
                if pr_table_X(k,1)~=0&&pr_table_X(m,1)~=0&&pr_table_combine(k,m)~=0
                    mi_sample=mi_sample+pr_table_combine(k,m)*log10(pr_table_combine(k,m)/(pr_table_X(k,1)*pr_table_Y(m,1)));
                end               
            end
        end
        % MI ������� ��� ǥ�� �Է��ϰ�, �ٽ� ����� ǥ���� �ʱ�ȭ
        mi_result(i,j)=mi_sample;
        pr_table_Y=zeros(num_bin,1);
        pr_table_combine=zeros(num_bin,num_bin);
        mi_sample=0;
        min_var_Y=0;
        max_var_Y=0;
        dif_var_Y=0;
        dif_cl_Y=0;
    end
    pr_table_X=zeros(num_bin,1);
    pr_tag=zeros(num_ano_data,1);
    min_var_X=0;
    max_var_X=0;
    dif_var_X=0;
    dif_cl_X=0;
end
% ���ʿ��� ������ �����մϴ�.
clear pr_table_X pr_tag pr_table_Y pr_table_combine min_var_X min_var_Y max_var_X max_var_Y dif_var_X dif_var_Y dif_cl_X dif_cl_Y mi_sample i j k l m n 
% Part2. ��.


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
    min_var_Y=min(ano_data_set(:,i));
    max_var_Y=max(ano_data_set(:,i));
    dif_var_Y=max_var_Y-min_var_Y;
    dif_cl_Y=dif_var_Y/num_bin;
    % Ym Ư�� �м�
    min_var_Ym=min(Ym(:,1));
    max_var_Ym=max(Ym(:,1));
    dif_var_Ym=max_var_Ym-min_var_Ym;
    dif_cl_Ym=dif_var_Ym/num_bin;
    % Ym �� Ȯ������ �ϼ�
    for j=1:num_ano_data-1
        n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
        % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
        pr_tag_Y(j,2)=n;
    end
    pr_t_Ym=pr_t_Ym/(num_ano_data-1);
    % Y�� Ym�� ����Ȯ������ �ϼ�
    for j=2:num_ano_data
        n=floor((ano_data_set(j,i)-min_var_Y)/dif_cl_Y)+1;
        if n==num_bin+1
            n=n-1;
        end
        % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
        pr_tag_Y(j-1,1)=n;
        pr_t_comb_Y(n,pr_tag_Y(j-1,2))=pr_t_comb_Y(n,pr_tag_Y(j-1,2))+1;
    end
    % ������ ���� ������� 0���� ũ�� 1���� ���� Ȯ�� ���� �˴ϴ�.
    pr_t_comb_Y=pr_t_comb_Y/(num_ano_data-1);
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
            min_var_X=min(Xt(:,1));
            max_var_X=max(Xt(:,1));
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            %XY����Ȯ�� �� XYYm����Ȯ�� �ϼ���Ű��
            for j=1:num_ano_data-tau
                n=floor((ano_data_set(j,m)-min_var_X)/dif_cl_X)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))=pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))+1;
                pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))=pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))+1;                
            end
            pr_t_comb_XY=pr_t_comb_XY/(num_ano_data-tau);
            pr_t_tri=pr_t_tri/(num_ano_data-tau);
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
clear i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym Xt
% Part3. ��.
% ��.


% Part4. 
% Monte Carlo Simulator �κ� 
% randsample(); �Լ��� �̿��Ͽ�, 
% mutual information�� transfer entropy�� ���� �ݺ��Ͽ� �����ϰ�,
% ������ ���ڱ��� ����� �����Ͽ�,(�⺻���� 1000ȸ, f=1000;)
% �� ���� ���� mean �� standard deviation�� ����մϴ�.
% �����  �Ʒ��� ���� ������İ� ���� ���·� ��Ÿ���ϴ�.
mi_mean_rand=zeros(num_var,num_var); %��ü ���
mi_sigma_rand=zeros(num_var,num_var); %��ü ǥ������
mi_sample_rand=zeros(1000,1); %������ mi ���
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
    for j=1:num_var
        Y=ano_data_set(:,j);
        for f=1:1000 %f�� ����Ƚ��
            rand_X=randsample(X,num_ano_data);        
            min_var_X=min(rand_X(:,1));
            max_var_X=max(rand_X(:,1));
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            % ���� X�� Ȯ������ǥ�� �ϼ��ϰ�, ������ ������ ���� tag �ϼ�
            % tag�� �ڿ� ����Ȯ�������� ����� �� ���˴ϴ�.
            for m=1:num_ano_data
                n=floor((rand_X(m,1)-min_var_X)/dif_cl_X)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_tag(m,1)=n;
                pr_table_X(n,1)=pr_table_X(n,1)+1;
            end
            pr_table_X=pr_table_X/num_ano_data;
            % X�� Ȯ�������� �ϼ��Ǹ�, ������ sink ����(Y)������ MI�� ��� ����ϴ� ����
            rand_Y=randsample(Y,num_ano_data);
            min_var_Y=min(rand_Y(:,1));
            max_var_Y=max(rand_Y(:,1));
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Y�� Ȯ�������� X,Y ����Ȯ������ �ϼ�
            for m=1:num_ano_data
                n=floor((rand_Y(m,1)-min_var_Y)/dif_cl_Y)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_table_Y(n,1)=pr_table_Y(n,1)+1;
                pr_table_combine(pr_tag(m,1),n)=pr_table_combine(pr_tag(m,1),n)+1;
            end
            pr_table_Y=pr_table_Y/num_ano_data;
            pr_table_combine=pr_table_combine/num_ano_data;
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
        mi_sample_rand=zeros(1000,1);
    end
end
% ���ʿ��� ������ �����մϴ�.
clear pr_table_X pr_tag pr_table_Y pr_table_combine min_var_X min_var_Y max_var_X max_var_Y dif_var_X dif_var_Y dif_cl_X dif_cl_Y mi_sample i j k l m n X Y rand_X rand_Y mi_sample_rand mi_result_rand
% Part2. �̿��Ͽ� mutual information�� Monte Carlo Simulation ��.
% Part3. �̿��Ͽ� transfer entropy�� Monte Carlo Simulation ����.
% ����� �Ʒ��� ���� ��������� ���·� ��Ÿ���ϴ�.
te_mean_rand=zeros(num_var,num_var,max_tau);
te_sigma_rand=zeros(num_var,num_var,max_tau);
te_sample_rand=zeros(1000,max_tau);
% te ��꿡 ���Ǵ� �����
te_sample=0;
% Sink ���(Y)�� i�� �����Ͽ� ����
for i=1:num_var
    % source ���(X)�� m���� �����Ͽ� ����
    for m=1:num_var
        X=ano_data_set(:,m);
        Y=ano_data_set(:,i);
        for f=1:1000
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
            min_var_Y=min(rand_Y(:,1));
            max_var_Y=max(rand_Y(:,1));
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Ym Ư�� �м�
            min_var_Ym=min(Ym(:,1));
            max_var_Ym=max(Ym(:,1));
            dif_var_Ym=max_var_Y-min_var_Y;
            dif_cl_Ym=dif_var_Y/num_bin;
            % Ym �� Ȯ������ �ϼ�
            for j=1:num_ano_data-1
                n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
                % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
                pr_tag_Y(j,2)=n;
            end
            pr_t_Ym=pr_t_Ym/(num_ano_data-1);
            % Y�� Ym�� ����Ȯ������ �ϼ�
            for j=1:num_ano_data-1
                n=floor((rand_Y(j+1,1)-min_var_Y)/dif_cl_Y)+1;
                if n==num_bin+1
                    n=n-1;
                end
                % tag�� ���ܼ�, ���߿� ����Ȯ������ ��� ��, �ٽ� ���
                pr_tag_Y(j,1)=n;
                pr_t_comb_Y(n,pr_tag_Y(j,2))=pr_t_comb_Y(n,pr_tag_Y(j,2))+1;
            end
            % ������ ���� ������� 0���� ũ�� 1���� ���� Ȯ�� ���� �˴ϴ�.
            pr_t_comb_Y=pr_t_comb_Y/(num_ano_data-1);
            % �����ϴ� �ð��Ը� tau�� 1���� max_tau���� ����
            for tau=1:max_tau
                % tau ��ŭ�� �ð����� ���� sink���(X) ~ m���� �� Ym�� ����Ȯ������ �ϼ�
                % X���� Ư�� ����
                Xt=zeros(num_ano_data-tau,1);
                for j=1:num_ano_data-tau
                    Xt(j,1)=ano_data_set(j,m);
                end                
                min_var_X=min(Xt(:,1));
                max_var_X=max(Xt(:,1));
                dif_var_X=max_var_X-min_var_X;
                dif_cl_X=dif_var_X/num_bin;
                for j=1:num_ano_data-tau
                    n=floor((rand_X(j,1)-min_var_X)/dif_cl_X)+1;
                    if n==num_bin+1
                        n=n-1;
                    end
                    pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))=pr_t_comb_XY(n,pr_tag_Y(j+tau-1,2))+1;
                    pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))=pr_t_tri(n, pr_tag_Y(j+tau-1,1), pr_tag_Y(j+tau-1,2))+1;                
                end
                pr_t_comb_XY=pr_t_comb_XY/(num_ano_data-tau);
                pr_t_tri=pr_t_tri/(num_ano_data-tau);
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
        te_sample_rand=zeros(1000,max_tau);
    end
end
% ���ʿ��� ������ �����մϴ�.
clear f i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym X Y rand_X rand_Y te_sample_rand num_data_day num_day timestep period_avg Xt
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
mi_threshold=0;
te_threshold=0;
% X(source)�� i, Y(sink)�� j�� rotation
for i=1:num_var
    for j=1:num_var
        mi_threshold=mi_mean_rand(i,j)+c*mi_sigma_rand(i,j);
        if mi_final(i,j) < mi_threshold
            mi_final(i,j)=0;
        end
        for k=1:max_tau
            te_threshold=te_mean_rand(i,j,k)+c*te_sigma_rand(i,j);
            if te_final(i,j,k) < te_threshold
                te_final(i,j,k)=0;
            end
        end
    end
end
clear c i j k mi_threshold te_threshold mi_mean_rand mi_sigma_rand te_mean_rand te_sigma_rand
% Part5. ��.


% Part6.
% normalizing �κ� : mi�� te�� max ���� log10(num_bin)���� ��ü�� ������
% ����� mi_final_nom�� te_final_nom ����
mi_final_nom=mi_final/num_bin;
te_final_nom=te_final/num_bin;
% Part6. ��.


% Part7.
% �� �������� �ִ밪�� ���ܼ� ������� ����� �κ�
% ����� adj_magnitude(num_var,num_var)�� adj_time(num_var,num_var) ����
% ����� ���� ����
adj_magnitude=zeros(num_var,num_var);
adj_time=zeros(num_var,num_var);
teadd=0;
for i=1:num_var
    for j=1:num_var
        teadd=find(te_final_nom(i,j,:)==max(te_final_nom(i,j,:)));
        if mi_final_nom(i,j)<te_final_nom(i,j,teadd)
            adj_magnitude(i,j)=te_final_nom(i,j,teadd);
            adj_time(i,j)=teadd;
        else
            adj_magnitude(i,j)=mi_final_nom(i,j);
            adj_time(i,j)=0;
        end
    end
end
clear i j k teadd 
% Part7. ��.

% ���� ����� ���̻� �ʿ���ٸ� �Ʒ��κб��� �����ϸ�, ����� �����ϴ�.
clear num_data_day num_day timestep period_avg
clear num_var num_ano_data num_bin max_tau