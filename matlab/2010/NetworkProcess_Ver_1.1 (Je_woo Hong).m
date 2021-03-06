% 이 프로그램의 작성은 현재진행중 입니다. (made by Je-woo Hong)
% 의견 주시면 반영하여 열심히 만들겠습니다.

% Version 1.0 (2011-03-09, Je-woo Hong)
% 현재는 data를 matlab 내부에서 정의하여 옮겨와서 사용해야 합니다.
% 결과 값들도, matlab 내부에 출력됩니다.
% 앞으로의 업그레이드 방향
% - data 및 result를 엑셀(또는 notepad)에서 읽고 쓰기
% - 자동으로 한달씩 끊어서 계산하기
% - 러델&쿠마 2번째 논문의 변수들 계산
% - 러델&쿠마 1번째 논문의 몬테카를로 시뮬레이션 및 인접행렬 등 추가연산 추가(Ver1.1에서 해결)

% Version 1.1 (2011-04-28, Je-woo Hong)
% 추가 및 수정 된 사항
% - Part4. MonteCarlo Simulation 추가
% - Part3. Transfer Entropy 계산 시 오류(data 시간순서) 수정
% - Part5. Weighted Cut 추가
% - Part6. Normalizing 추가
% - Part7. 인접행렬 추가


% Part1.
% data_set을 anomaly로 만드는 부분입니다.
% 가장먼저 data_set을 정의하고, data를 입력하세요.
% ex) data_set=zeros(데이터수, 변수수);
% 먼저 anomaly를 만들기 위해 평균하는 기간을 지정합니다. "period_avg"
% 기본값은 5일 (5일(오늘과 지난 4일)동안 같은시간의 데이터를 평균내어 사용)
period_avg=5;
% 데이터의 "timestep" 결정 (단위 : sec), 기본값은 30min = 1800sec
timestep=1800;
% 데이터의 기본 정보 값을 조사합니다.
size_data=size(data_set);
num_var=size_data(1,2)
num_data=size_data(1,1)
num_data_day=24*3600/timestep
num_day=num_data(1,1)/num_data_day
avg_sample=0;
ano_data_set=zeros(num_data-num_data_day*(period_avg-1),num_var);
% anomaly 데이터를 형성합니다.
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
% anomaly 데이터 정보를 조사합니다.
size_data=size(ano_data_set);
num_ano_data=size_data(1,1)
num_day=num_ano_data(1,1)/num_data_day
% 뒤에 불필요한 변수를 제거합니다.
clear i j k avg_sample size_data num_data
% Part1. 끝.


% Part2.
% anomaly 데이터를 이용하여, 
% 변수끼리 가질 수 있는 모든 경우의, mutual information을 계산합니다.
% 결과는 mi_result 라는 이름으로 인접행렬의 형태로 나타납니다.
% 먼저 bin 수를 "num_bin"에 정의해주세요. 기본값은 11 (러델&쿠마)
num_bin=11;
% 확률분포 계산 및 MI계산에 필요한 변수 정의
pr_table_X=zeros(num_bin,1);
pr_table_Y=zeros(num_bin,1);
pr_table_combine=zeros(num_bin,num_bin);
pr_tag=zeros(num_ano_data,1);
mi_result=zeros(num_var,num_var);
mi_sample=0;
% Source 변수(X)를 i로 정의하여 루프
for i=1:num_var 
    min_var_X=min(ano_data_set(:,i));
    max_var_X=max(ano_data_set(:,i));
    dif_var_X=max_var_X-min_var_X;
    dif_cl_X=dif_var_X/num_bin;
    % 먼저 X의 확률분포표를 완성하고, 데이터 순번에 따른 tag 완성
    % tag는 뒤에 결합확률분포를 계산할 때 사용됩니다.
    for m=1:num_ano_data
        n=floor((ano_data_set(m,i)-min_var_X)/dif_cl_X)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_tag(m,1)=n;
        pr_table_X(n,1)=pr_table_X(n,1)+1;
    end
    pr_table_X=pr_table_X/num_ano_data;
    % X의 확률분포가 완성되면, 나머지 sink 변수(Y)에대한 MI를 모두 계산하는 루프
    for j=1:num_var
        min_var_Y=min(ano_data_set(:,j));
        max_var_Y=max(ano_data_set(:,j));
        dif_var_Y=max_var_Y-min_var_Y;
        dif_cl_Y=dif_var_Y/num_bin;
        % Y의 확률분포와 X,Y 결합확률분포 완성
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
        % MI 계산 루프
        for k=1:num_bin
            for m=1:num_bin
                if pr_table_X(k,1)~=0&&pr_table_X(m,1)~=0&&pr_table_combine(k,m)~=0
                    mi_sample=mi_sample+pr_table_combine(k,m)*log10(pr_table_combine(k,m)/(pr_table_X(k,1)*pr_table_Y(m,1)));
                end               
            end
        end
        % MI 계산결과를 결과 표에 입력하고, 다시 사용할 표들은 초기화
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
% 불필요한 변수를 제거합니다.
clear pr_table_X pr_tag pr_table_Y pr_table_combine min_var_X min_var_Y max_var_X max_var_Y dif_var_X dif_var_Y dif_cl_X dif_cl_Y mi_sample i j k l m n 
% Part2. 끝.


% Part3.
% anomaly 데이터를 이용하여, 
% 지정한 시간규모(max_tau) 까지의,
% 변수끼리 가질 수 있는 모든 경우에 대한, transfer entropy를 계산합니다.
% 결과는 te_result 라는 이름으로 인접행렬의 형태로 나타납니다.
% 먼저 조사하고자 하는 시간규모 한계인 "max_tau"를 지정합니다. 
% 기본값 36 (30분 간격일 때 18시간, 러델&쿠마)
max_tau=36;
% te 계산에 사용되는 빈공간
te_result=zeros(num_var,num_var,max_tau);
te_sample=0;
% Sink 노드(Y)를 i로 정의하여 루프
for i=1:num_var
    pr_t_comb_Y=zeros(num_bin,num_bin); % 변수 순서는 Y, Ym
    % Ym 의 의미는 Y의 1개 timestep 전의 데이터열을 뜻합니다.(Y_t-dt)
    % 계산에 필요한 표들을 초기화
    pr_t_Ym=zeros(num_bin,1);
    pr_t_comb_XY=zeros(num_bin,num_bin); %X, Y
    pr_t_tri=zeros(num_bin,num_bin,num_bin); %X, Y, Ym
    pr_tag_Y=zeros(num_ano_data-1,2); %Y, Ym
    % Y와 -1 timestep 차이가 나는 데이터(Ym) 형성
    % Ym의 첫번째 데이터는 Ym의 t=2일 때 데이터부터 시작됨.(t=1은 없음)
    Ym=zeros(num_ano_data-1,1);
    for j=1:num_ano_data-1
        Ym(j,1)=ano_data_set(j,i);
    end
    % Y 변수 특성 분석
    min_var_Y=min(ano_data_set(:,i));
    max_var_Y=max(ano_data_set(:,i));
    dif_var_Y=max_var_Y-min_var_Y;
    dif_cl_Y=dif_var_Y/num_bin;
    % Ym 특성 분석
    min_var_Ym=min(Ym(:,1));
    max_var_Ym=max(Ym(:,1));
    dif_var_Ym=max_var_Ym-min_var_Ym;
    dif_cl_Ym=dif_var_Ym/num_bin;
    % Ym 의 확률분포 완성
    for j=1:num_ano_data-1
        n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
        if n==num_bin+1
            n=n-1;
        end
        pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
        % tag를 남겨서, 나중에 결합확률분포 계산 때, 다시 사용
        pr_tag_Y(j,2)=n;
    end
    pr_t_Ym=pr_t_Ym/(num_ano_data-1);
    % Y와 Ym의 결합확률분포 완성
    for j=2:num_ano_data
        n=floor((ano_data_set(j,i)-min_var_Y)/dif_cl_Y)+1;
        if n==num_bin+1
            n=n-1;
        end
        % tag를 남겨서, 나중에 결합확률분포 계산 때, 다시 사용
        pr_tag_Y(j-1,1)=n;
        pr_t_comb_Y(n,pr_tag_Y(j-1,2))=pr_t_comb_Y(n,pr_tag_Y(j-1,2))+1;
    end
    % 데이터 수로 나눠줘야 0보다 크고 1보다 작은 확률 값이 됩니다.
    pr_t_comb_Y=pr_t_comb_Y/(num_ano_data-1);
    % tau 만큼의 시간차를 갖는 sink노드(X) ~ m루프 와 Ym의 결합확률분포 완성
    % source 노드(X)를 m으로 정의하여 루프
    for m=1:num_var
        % 조사하는 시간규모 tau를 1부터 max_tau까지 조사
        for tau=1:max_tau
            % X변수 특성 조사
            Xt=zeros(num_ano_data-tau,1);
            for j=1:num_ano_data-tau
                Xt(j,1)=ano_data_set(j,m);
            end                
            min_var_X=min(Xt(:,1));
            max_var_X=max(Xt(:,1));
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            %XY결합확률 및 XYYm결합확률 완성시키기
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
            % te 계산 루프
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
            % 결과를 입력하고, 계산에 반복 사용되는 표 초기화
            te_result(m,i,tau)=te_sample;
            pr_t_tri=zeros(num_bin,num_bin,num_bin);
            pr_t_comb_XY=zeros(num_bin,num_bin);
            te_sample=0;     
        end
    end
end
% 불필요한 변수를 제거합니다.
clear i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym Xt
% Part3. 끝.
% 끝.


% Part4. 
% Monte Carlo Simulator 부분 
% randsample(); 함수를 이용하여, 
% mutual information과 transfer entropy의 값을 반복하여 추출하고,
% 지정한 숫자까지 계산을 수행하여,(기본값은 1000회, f=1000;)
% 각 값에 대한 mean 과 standard deviation을 기록합니다.
% 결과는  아래와 같이 인접행렬과 같은 형태로 나타납니다.
mi_mean_rand=zeros(num_var,num_var); %전체 평균
mi_sigma_rand=zeros(num_var,num_var); %전체 표준편차
mi_sample_rand=zeros(1000,1); %샘플한 mi 목록
% Part2. 이용해서, mutual information, Monte Carlo simulation 수행
% 확률분포 계산 및 MI계산에 필요한 변수 정의
pr_table_X=zeros(num_bin,1);
pr_table_Y=zeros(num_bin,1);
pr_table_combine=zeros(num_bin,num_bin);
pr_tag=zeros(num_ano_data,1);
mi_sample=0;
X=zeros(num_ano_data,1);
Y=zeros(num_ano_data,1);
% Source 변수(X)를 i로 정의하여 루프
for i=1:num_var 
    X=ano_data_set(:,i);
    for j=1:num_var
        Y=ano_data_set(:,j);
        for f=1:1000 %f는 수행횟수
            rand_X=randsample(X,num_ano_data);        
            min_var_X=min(rand_X(:,1));
            max_var_X=max(rand_X(:,1));
            dif_var_X=max_var_X-min_var_X;
            dif_cl_X=dif_var_X/num_bin;
            % 먼저 X의 확률분포표를 완성하고, 데이터 순번에 따른 tag 완성
            % tag는 뒤에 결합확률분포를 계산할 때 사용됩니다.
            for m=1:num_ano_data
                n=floor((rand_X(m,1)-min_var_X)/dif_cl_X)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_tag(m,1)=n;
                pr_table_X(n,1)=pr_table_X(n,1)+1;
            end
            pr_table_X=pr_table_X/num_ano_data;
            % X의 확률분포가 완성되면, 나머지 sink 변수(Y)에대한 MI를 모두 계산하는 루프
            rand_Y=randsample(Y,num_ano_data);
            min_var_Y=min(rand_Y(:,1));
            max_var_Y=max(rand_Y(:,1));
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Y의 확률분포와 X,Y 결합확률분포 완성
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
            % MI 계산 루프
            for k=1:num_bin
                for m=1:num_bin
                    if pr_table_X(k,1)~=0&&pr_table_X(m,1)~=0&&pr_table_combine(k,m)~=0
                        mi_sample=mi_sample+pr_table_combine(k,m)*log10(pr_table_combine(k,m)/(pr_table_X(k,1)*pr_table_Y(m,1)));
                    end               
                end
            end
            % MI 계산결과를 결과 표에 입력하고, 다시 사용할 표들은 초기화
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
% 불필요한 변수를 제거합니다.
clear pr_table_X pr_tag pr_table_Y pr_table_combine min_var_X min_var_Y max_var_X max_var_Y dif_var_X dif_var_Y dif_cl_X dif_cl_Y mi_sample i j k l m n X Y rand_X rand_Y mi_sample_rand mi_result_rand
% Part2. 이용하여 mutual information의 Monte Carlo Simulation 끝.
% Part3. 이용하여 transfer entropy의 Monte Carlo Simulation 시작.
% 결과는 아래와 같은 인접행렬의 형태로 나타납니다.
te_mean_rand=zeros(num_var,num_var,max_tau);
te_sigma_rand=zeros(num_var,num_var,max_tau);
te_sample_rand=zeros(1000,max_tau);
% te 계산에 사용되는 빈공간
te_sample=0;
% Sink 노드(Y)를 i로 정의하여 루프
for i=1:num_var
    % source 노드(X)를 m으로 정의하여 루프
    for m=1:num_var
        X=ano_data_set(:,m);
        Y=ano_data_set(:,i);
        for f=1:1000
            rand_X=randsample(X,num_ano_data);
            rand_Y=randsample(Y,num_ano_data);
            % Y와 1 timestep 차이가 나는 데이터(Ym) 형성
            % Ym의 의미는 Y의 1개 timestep 전의 데이터열을 뜻합니다.(Y_t-dt) Ym은 t=2부터 존재
            Ym=zeros(num_ano_data-1,1);
            for j=1:num_ano_data-1
                Ym(j,1)=rand_Y(j,1);
            end
            % 계산에 필요한 표들을 초기화
            pr_t_comb_Y=zeros(num_bin,num_bin); % 변수 순서는 Y, Ym
            pr_t_Ym=zeros(num_bin,1);
            pr_t_comb_XY=zeros(num_bin,num_bin); %X, Y
            pr_t_tri=zeros(num_bin,num_bin,num_bin); %X, Y, Ym
            pr_tag_Y=zeros(num_ano_data-1,2); %Y, Ym            
            % Y 변수 특성 분석
            min_var_Y=min(rand_Y(:,1));
            max_var_Y=max(rand_Y(:,1));
            dif_var_Y=max_var_Y-min_var_Y;
            dif_cl_Y=dif_var_Y/num_bin;
            % Ym 특성 분석
            min_var_Ym=min(Ym(:,1));
            max_var_Ym=max(Ym(:,1));
            dif_var_Ym=max_var_Y-min_var_Y;
            dif_cl_Ym=dif_var_Y/num_bin;
            % Ym 의 확률분포 완성
            for j=1:num_ano_data-1
                n=floor((Ym(j,1)-min_var_Ym)/dif_cl_Ym)+1;
                if n==num_bin+1
                    n=n-1;
                end
                pr_t_Ym(n,1)=pr_t_Ym(n,1)+1;
                % tag를 남겨서, 나중에 결합확률분포 계산 때, 다시 사용
                pr_tag_Y(j,2)=n;
            end
            pr_t_Ym=pr_t_Ym/(num_ano_data-1);
            % Y와 Ym의 결합확률분포 완성
            for j=1:num_ano_data-1
                n=floor((rand_Y(j+1,1)-min_var_Y)/dif_cl_Y)+1;
                if n==num_bin+1
                    n=n-1;
                end
                % tag를 남겨서, 나중에 결합확률분포 계산 때, 다시 사용
                pr_tag_Y(j,1)=n;
                pr_t_comb_Y(n,pr_tag_Y(j,2))=pr_t_comb_Y(n,pr_tag_Y(j,2))+1;
            end
            % 데이터 수로 나눠줘야 0보다 크고 1보다 작은 확률 값이 됩니다.
            pr_t_comb_Y=pr_t_comb_Y/(num_ano_data-1);
            % 조사하는 시간규모 tau를 1부터 max_tau까지 조사
            for tau=1:max_tau
                % tau 만큼의 시간차를 갖는 sink노드(X) ~ m루프 와 Ym의 결합확률분포 완성
                % X변수 특성 조사
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
                % te 계산 루프
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
                % 결과를 입력하고, 계산에 반복 사용되는 표 초기화
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
% 불필요한 변수를 제거합니다.
clear f i j m n p s pr_t_Ym pr_t_comb_XY pr_t_comb_Y pr_t_tri te_sample tau dif_cl_X dif_cl_Y dif_var_X dif_var_Y min_var_X min_var_Y max_var_X max_var_Y pr_tag_Y min_var_Ym max_var_Ym dif_var_Ym dif_cl_Ym Ym X Y rand_X rand_Y te_sample_rand num_data_day num_day timestep period_avg Xt
% Part3. 이용 transfer entropy의 Monte Carlo Simulation 끝.
% Part4. 끝.


% Part5.
% Weighted Cut 부분 : 통계적 유의성(mu+sigma)과 비교하여, 작을 경우 제거
%   T > "mu(T_ss)+c*sigma(T_ss) = threshold of T"
%   95% : c=1.66
%   99% : c=2.36
% 여기서는 기본값 c=1.66 (95%) 사용
c=1.66;
% 결과로 mi_final(num_var,num_var) 과 te_final(num_var,num_var,max_tau) 을 생산
mi_final=mi_result;
te_final=te_result;
% 사용되는 빈 변수 지정
mi_threshold=0;
te_threshold=0;
% X(source)는 i, Y(sink)는 j로 rotation
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
% Part5. 끝.


% Part6.
% normalizing 부분 : mi와 te의 max 값인 log10(num_bin)으로 전체를 나눠줌
% 결과로 mi_final_nom과 te_final_nom 생산
mi_final_nom=mi_final/num_bin;
te_final_nom=te_final/num_bin;
% Part6. 끝.


% Part7.
% 각 변수간의 최대값만 남겨서 인접행렬 만드는 부분
% 결과로 adj_magnitude(num_var,num_var)와 adj_time(num_var,num_var) 생산
% 사용할 변수 지정
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
% Part7. 끝.

% 만약 계산이 더이상 필요없다면 아래부분까지 실행하면, 결과만 남습니다.
clear num_data_day num_day timestep period_avg
clear num_var num_ano_data num_bin max_tau