    %% Ogive test for "kinematic sensible heat flux"
    % History
    % 2012-11, J-W Hong, for Je-Ju
    % 2013-06-01, J-W Hong, for Hong-Chun and Sam-Chuk, using 1day TOA5 data


    %% Information
    length = 24; % data length in hour
    Hz = 10; % sampling rate
    po_u = 2; % column number of u
    avg_period=[1;2;3;4;5;6;8;10;12;15;20;24;30;40;60;120]; 
        % averaging period in [min] for Ogive-test each 2 hours
    temp = size(avg_period);
    n_avg = temp(1,1);
    Ogive=zeros(n_avg,length/2); 
    clear temp


    %% Sprit 1 day data into 2 hours

    sample = zeros(length/2,Hz*2*3600,4);
    for i = 1:length/2
        for j = 1:Hz*2*3600
            for k = 1:4
                sample(i,j,k) = data(j+(i-1)*Hz*2*3600,po_u -1 +k);
            end
        end
    end
    clear i j k


    %% Ogive Test

    for i = 1:length/2
        for j = 1:n_avg
            n_data = avg_period(j,1)*60*Hz; 
            % number of raw data for single flux point
            n_point = 120/avg_period(j,1); 
            temp_flux=zeros(n_point,1); % bin for result
            for k = 1:n_point
                temp_cov=0;	
                temp_u=zeros(n_data,1);
                temp_v=zeros(n_data,1);
                temp_w=zeros(n_data,1);
                temp_ts=zeros(n_data,1);
                % extract for single flux
                for l = 1:n_data 
                    temp_u(l,1)=sample(i,n_data*(k-1)+l,1);
                    temp_v(l,1)=sample(i,n_data*(k-1)+l,2);
                    temp_w(l,1)=sample(i,n_data*(k-1)+l,3);
                    temp_ts(l,1)=sample(i,n_data*(k-1)+l,4);
                end
                u_bar = mean(temp_u); % for double rotation
                v_bar = mean(temp_v);
                w_bar = mean(temp_w);
                alpha = atan(v_bar/u_bar);
                beta = atan(w_bar/(sqrt(u_bar^2+v_bar^2)));
                for l = 1:n_data
                    temp_u(l,1) = cos(beta)*(cos(alpha)*temp_u(l,1)+sin(alpha)*temp_v(l,1))+sin(beta)*temp_w(l,1);
                    temp_v(l,1) = -sin(alpha)*temp_u(l,1)+cos(alpha)*temp_v(l,1);
                    temp_w(l,1) = -sin(beta)*(cos(alpha)*temp_u(l,1)+sin(alpha)*temp_v(l,1))+cos(beta)*temp_w(l,1);
                end
                w_mean = mean(temp_w);
                ts_mean = mean(temp_ts);
                for l = 1:n_data
                    temp_w(l,1)=temp_w(l,1)-w_mean;
                    temp_ts(l,1)=temp_ts(l,1)-ts_mean;
                end
                for l = 1:n_data
                    temp_cov=temp_cov + temp_w(l,1)*temp_ts(l,1);
                end
                temp_flux(k,1)=temp_cov/n_data;
            end
            Ogive(j,i)=mean(temp_flux);
        end
    end

    clear w_mean ts_mean n_data n_point temp_flux i j k l Hz alpha avg_period beta length n_avg po_u sample temp_cov temp_ts temp_u temp_v temp_w u_bar v_bar w_bar
    clear textdata
