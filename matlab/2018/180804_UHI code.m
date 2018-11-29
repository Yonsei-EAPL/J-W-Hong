%% UHI SE and KP sites
% 1-hourly from 1983-01-01 01:00 ~ 2018-01-01 00:00

%% data x
% 1: date (excel num)
% 2: hour
% 3: T KP site
% 4: T SE site
% 5: UHIi : T SE-KP
% 6: precipitation (mm/hour)
% 7: month

%% data of a day
l = 306816; % total length
temp1 = 0; %n_days
temp2 = 0; %temp
for i = 1:l
    if i==1
        temp1=1;
        temp2=x(i,1);
    else
        if temp2 ~=x(i,1);
            temp1 = temp1+1; %n_days
            temp2 = x(i,1);
        end
    end
end
clear i temp2

y=zeros(temp1,20); % output
m=temp1; %n_days
clear temp1
% 1: date (excel num)
% 2: n of a day KP site
% 3: n of a day SE site
% 4: Tmax KP
% 5: Tmax hour KP
% 6: Tmin KP
% 7: Tmin hour KP
% 8: Tmean KP
% 9: DTR(daily temperature range) KP
% 10: Tmax SE
% 11: Tmax hour SE
% 12: Tmin SE
% 13: Tmin hour SE
% 14: Tmean SE
% 15: DTR(daily temperature range) SE
% 16: UHIi min
% 17: UHIi min hour
% 18: UHIi max
% 19: UHIi max hour
% 20: precipitation (mm/day)

%% 
temp1 = 0; %n_days
temp2 = 0; %temp
for i = 1:l
    if i==1 % first data
        temp1=1;
        temp2=x(i,1);
        y(temp1,1) = temp2;
        % KP
        if x(i,3)>-90
            y(temp1,2) = y(temp1,2)+1;
            y(temp1,4) = x(i,3);
            y(temp1,5) = x(i,2);
            y(temp1,6) = x(i,3);
            y(temp1,7) = x(i,2);
            y(temp1,8) = y(temp1,8)+x(i,3);
        end
        % SE
        if x(i,4)>-90
            y(temp1,3) = y(temp1,3)+1;
            y(temp1,10) = x(i,4);
            y(temp1,11) = x(i,2);
            y(temp1,12) = x(i,4);
            y(temp1,13) = x(i,2);
            y(temp1,14) = y(temp1,14)+ x(i,4);
        end
        % UHIi
        if x(i,5)>-90
            y(temp1,16) = x(i,5);
            y(temp1,17) = x(i,2);
            y(temp1,18) = x(i,5);
            y(temp1,19) = x(i,2);
        end
        % precipitation
        if x(i,6)>0
            y(temp1,20) = y(temp1,20)+ x(i,6);
        end
    else
        if temp2 ~=x(i,1); % day change
            temp1 = temp1+1; %n_days
            temp2 = x(i,1);
            y(temp1,1) = temp2;
            % statistics
            % KP
            if y(temp1-1,2)>1
                y(temp1-1,8)= y(temp1-1,8)/y(temp1-1,2); % mean
                y(temp1-1,9)=y(temp1-1,4)-y(temp1-1,6); % DTR
            end
            % SE
            if y(temp1-1,3)>1
                y(temp1-1,14)= y(temp1-1,14)/y(temp1-1,3); % mean
                y(temp1-1,15)=y(temp1-1,10)-y(temp1-1,12); % DTR
            end
            
            % KP
            if x(i,3)>-90
                y(temp1,2) = y(temp1,2)+1;
                y(temp1,4) = x(i,3);
                y(temp1,5) = x(i,2);
                y(temp1,6) = x(i,3);
                y(temp1,7) = x(i,2);
                y(temp1,8) = y(temp1,8)+x(i,3);
            end
            % SE
            if x(i,4)>-90
                y(temp1,3) = y(temp1,3)+1;
                y(temp1,10) = x(i,4);
                y(temp1,11) = x(i,2);
                y(temp1,12) = x(i,4);
                y(temp1,13) = x(i,2);
                y(temp1,14) = y(temp1,14)+ x(i,4);
            end
            % UHIi
            if x(i,5)>-90
                y(temp1,16) = x(i,5);
                y(temp1,17) = x(i,2);
                y(temp1,18) = x(i,5);
                y(temp1,19) = x(i,2);
            end
            % precipitation
            if x(i,6)>0
                y(temp1,20) = y(temp1,20)+ x(i,6);
            end            
        else % normal
            % KP
            if x(i,3)>-90
                y(temp1,2) = y(temp1,2)+1;
                if x(i,3)>y(temp1,4)
                    y(temp1,4) = x(i,3);
                    y(temp1,5) = x(i,2);
                end
                if x(i,3)<y(temp1,6)
                    y(temp1,6) = x(i,3);
                    y(temp1,7) = x(i,2);
                end
                y(temp1,8) = y(temp1,8)+x(i,3);
            end
            % SE
            if x(i,4)>-90
                y(temp1,3) = y(temp1,3)+1;
                if x(i,4)>y(temp1,10)
                    y(temp1,10) = x(i,4);
                    y(temp1,11) = x(i,2);
                end
                if x(i,4)<y(temp1,12)
                    y(temp1,12) = x(i,4);
                    y(temp1,13) = x(i,2);
                end
                y(temp1,14) = y(temp1,14)+ x(i,4);
            end
            % UHIi
            if x(i,5)>-90
                if x(i,5)<y(temp1,16)
                    y(temp1,16) = x(i,5);
                    y(temp1,17) = x(i,2);
                end
                if x(i,5)>y(temp1,18)
                    y(temp1,18) = x(i,5);
                    y(temp1,19) = x(i,2);
                end
            end
            % precipitation
            if x(i,6)>0
                y(temp1,20) = y(temp1,20)+ x(i,6);
            end
        end
    end
end
clear i


