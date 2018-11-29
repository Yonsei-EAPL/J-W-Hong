% data_UHI
% 1: YYYYMMDD00
% 2: site number
% 3: daily maximum UHI (T_AWS - T_±èÆ÷°øÇ×)
% 4: LCZ at maximum UHI
% 5: HH (time)
% 6: number of all data (24/23 are OK)
% 7: number of good data
% 8: number of gaps
% 9: rainy day(1/0)
% 10: hour after rainfall
% 11: solar radiation after rainfall
% 12: UHI final with environmental lapse rate

elr = 6.49/1000; % oC/m
kp = 11.38; % height of kimpo airport
data_uhi(1,12)=0;
for i = 1:length(data_uhi(:,1))
    % 110: 11.38m
    if data_uhi(i,2)==108 % 85.67
        data_uhi(i,12) = data_uhi(i,3)  + (85.67-kp)*elr;
    elseif data_uhi(i,2)==400 % 59.24
        data_uhi(i,12) = data_uhi(i,3)  + (59.24-kp)*elr;
    elseif data_uhi(i,2)==401 % 35.53 before 20160912/ 33.05 from 20160913
        if data_uhi(i,1)<2016091299
            data_uhi(i,12) = data_uhi(i,3)  + (35.53-kp)*elr;
        else
            data_uhi(i,12) = data_uhi(i,3)  + (33.05-kp)*elr;
        end
    elseif data_uhi(i,2)==402 % 56.21
        data_uhi(i,12) = data_uhi(i,3)  + (56.21-kp)*elr;
    elseif data_uhi(i,2)==403 % 58.26
        data_uhi(i,12) = data_uhi(i,3)  + (58.26-kp)*elr;
    elseif data_uhi(i,2)==404 % 80.40 before 20170928/ 13 from 20190929
        if data_uhi(i,1)<2017092899
            data_uhi(i,12) = data_uhi(i,3)  + (80.40-kp)*elr;
        else
            data_uhi(i,12) = data_uhi(i,3)  + (13-kp)*elr;
        end
    elseif data_uhi(i,2)==405 % 11.93
        data_uhi(i,12) = data_uhi(i,3)  + (11.93-kp)*elr;
    elseif data_uhi(i,2)==406 % 56.65
        data_uhi(i,12) = data_uhi(i,3)  + (56.65-kp)*elr;
    elseif data_uhi(i,2)==407 % 52.05
        data_uhi(i,12) = data_uhi(i,3)  + (52.05-kp)*elr;
    elseif data_uhi(i,2)==408 % 53.96
        data_uhi(i,12) = data_uhi(i,3)  + (53.96-kp)*elr;
    elseif data_uhi(i,2)==409 % 39.09
        data_uhi(i,12) = data_uhi(i,3)  + (39.09-kp)*elr;
    elseif data_uhi(i,2)==410 % 36.44
        data_uhi(i,12) = data_uhi(i,3)  + (36.44-kp)*elr;
    elseif data_uhi(i,2)==411 % 24.00 before 20100815/ 100.67 after 20100816
        if data_uhi(i,1)<2010081599
            data_uhi(i,12) = data_uhi(i,3)  + (24.00-kp)*elr;
        else
            data_uhi(i,12) = data_uhi(i,3)  + (100.67-kp)*elr;
        end
    elseif data_uhi(i,2)==412 % 103.08
        data_uhi(i,12) = data_uhi(i,3)  + (103.08-kp)*elr;
    elseif data_uhi(i,2)==413 % 28.63
        data_uhi(i,12) = data_uhi(i,3)  + (28.63-kp)*elr;
    elseif data_uhi(i,2)==414 % 128.62
        data_uhi(i,12) = data_uhi(i,3)  + (128.62-kp)*elr;
    elseif data_uhi(i,2)==415 % 31.73
        data_uhi(i,12) = data_uhi(i,3)  + (31.73-kp)*elr;
    elseif data_uhi(i,2)==416 % 57.92
        data_uhi(i,12) = data_uhi(i,3)  + (57.92-kp)*elr;
    elseif data_uhi(i,2)==417 % 45.00
        data_uhi(i,12) = data_uhi(i,3)  + (45-kp)*elr;
    elseif data_uhi(i,2)==418 % 10.66
        data_uhi(i,12) = data_uhi(i,3)  + (10.66-kp)*elr;
    elseif data_uhi(i,2)==419 % 266.44
        data_uhi(i,12) = data_uhi(i,3)  + (266.44-kp)*elr;
    elseif data_uhi(i,2)==421 % 34.73
        data_uhi(i,12) = data_uhi(i,3)  + (34.73-kp)*elr;
    elseif data_uhi(i,2)==422 % 333.80
        data_uhi(i,12) = data_uhi(i,3)  + (333.80-kp)*elr;
    elseif data_uhi(i,2)==423 % 56.08
        data_uhi(i,12) = data_uhi(i,3)  + (56.08-kp)*elr;
    elseif data_uhi(i,2)==424 % 57.29
        data_uhi(i,12) = data_uhi(i,3)  + (57.29-kp)*elr;
    elseif data_uhi(i,2)==425 % 87.86
        data_uhi(i,12) = data_uhi(i,3)  + (87.86-kp)*elr;
    elseif data_uhi(i,2)==509 % 92.36 before 20040904/ 141.64 after 20040905
        if data_uhi(i,1)<2004090499
            data_uhi(i,12) = data_uhi(i,3)  + (92.36-kp)*elr;
        else
            data_uhi(i,12) = data_uhi(i,3)  + (141.64-kp)*elr;
        end
    elseif data_uhi(i,2)==510 % 25.38
        data_uhi(i,12) = data_uhi(i,3)  + (25.38-kp)*elr;
    elseif data_uhi(i,2)==889 % 16.23
        data_uhi(i,12) = data_uhi(i,3)  + (16.23-kp)*elr;
    end
end
clear i
clear elr kp

save('data_uhi.mat','data_uhi');