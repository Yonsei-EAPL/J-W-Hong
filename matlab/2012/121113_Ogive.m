% Ogive test for "kinematic sensible heat flux"

avg_period=[1;2;3;4;5;6;8;10;12;15;20;24;30;40;60;120]; % averaging period[min] for Ogive-test each 120min
result=zeros(16,6); % 16 means number of avg_period, and 6 means 120min bins from 6am to 6pm

for m = 1:6

	for i = 1:16

		no_data=avg_period(i,1)*1200; % number of single-data for single-flux point
		no_point=120/avg_period(i,1); % number of flux results during 120min
		temp_flux=zeros(no_point,1);

		for j = 1:no_point

			w_mean=0;
			ts_mean=0; 
			temp_cov=0;	
            temp_u=zeros(no_data,1);
            temp_v=zeros(no_data,1);
            temp_w=zeros(no_data,1);
            temp_ts=zeros(no_data,1);

			for k = 1:no_data

				temp_u(k,1)=data(20*6*60*60+(m-1)*20*60*60*2+(j-1)*no_data+k,1);
                temp_v(k,1)=data(20*6*60*60+(m-1)*20*60*60*2+(j-1)*no_data+k,2);
                temp_w(k,1)=data(20*6*60*60+(m-1)*20*60*60*2+(j-1)*no_data+k,3);
				temp_ts(k,1)=data(20*6*60*60+(m-1)*20*60*60*2+(j-1)*no_data+k,4);

			end

            u_bar = mean(temp_u);
            v_bar = mean(temp_v);
            w_bar = mean(temp_w);
            alpha = atan(v_bar/u_bar);
            beta = atan(w_bar/(sqrt(u_bar^2+v_bar^2)));

            for k = 1:no_data

                temp_u(k,1) = cos(beta)*(cos(alpha)*temp_u(k,1)+sin(alpha)*temp_v(k,1))+sin(beta)*temp_w(k,1);
                temp_v(k,1) = -sin(alpha)*temp_u(k,1)+cos(alpha)*temp_v(k,1);
                temp_w(k,1) = -sin(beta)*(cos(alpha)*temp_u(k,1)+sin(alpha)*temp_v(k,1))+cos(beta)*temp_w(k,1);

            end

			w_mean = mean(temp_w(:,1));
			ts_mean = mean(temp_ts(:,1));
	
			for k = 1:no_data
		
				temp_w(k,1)=temp_w(k,1)-w_mean;
				temp_ts(k,1)=temp_ts(k,1)-ts_mean;

			end

			for k = 1:no_data
			
				temp_cov=temp_cov + temp_w(k,1)*temp_ts(k,1);

			end

			temp_flux(j,1)=temp_cov/no_data;

		end

		result(i,m)=mean(temp_flux(:,1));

	end

end

clear w_mean ts_mean no_data no_point temp_flux i j k m ans
