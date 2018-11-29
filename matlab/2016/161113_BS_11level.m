%% input
% data_ws1,2,3,4
% data_wd1,2,3,4

% R = input('=')

%% 

for jw = 1:11
    WD = data_wd1(:,jw);
    WD2 = data_wd2(:,jw);
    WD3 = data_wd3(:,jw);
    WD4 = data_wd4(:,jw);    

    data = data_ws1(:,jw);
    data2 = data_ws2(:,jw);
    data3 = data_ws3(:,jw);
    data4 = data_ws4(:,jw);

    %%
    temp_max = 8;
    temp_min = 0;

    WD_backup = WD;
    data_backup = data;
    a = size(WD_backup);
    sizeb = a(1,1);
    clear a
    for i = 1:sizeb
        if WD_backup(i,1)*(-1)+90<0
            WD_backup(i,1) = WD_backup(i,1)*(-1)+450;
        else
            WD_backup(i,1) = WD_backup(i,1)*(-1)+90;
        end
    end
    for i = 1:sizeb
        if data_backup(i,1) > temp_max
            data_backup(i,1) =temp_max;
        elseif data_backup(i,1) < temp_min
            data_backup(i,1) = temp_min;
        end
    end

    WD_backup2 = WD2;
    data_backup2 = data2;
    a = size(WD_backup2);
    sizeb = a(1,1);
    clear a
    for i = 1:sizeb
        if WD_backup2(i,1)*(-1)+90<0
            WD_backup2(i,1) = WD_backup2(i,1)*(-1)+450;
        else
            WD_backup2(i,1) = WD_backup2(i,1)*(-1)+90;
        end
    end
    for i = 1:sizeb
        if data_backup2(i,1) > temp_max
            data_backup2(i,1) =temp_max;
        elseif data_backup2(i,1) < temp_min
            data_backup2(i,1) =temp_min;
        end
    end

    WD_backup3 = WD3;
    data_backup3 = data3;
    a = size(WD_backup3);
    sizeb = a(1,1);
    clear a
    for i = 1:sizeb
        if WD_backup3(i,1)*(-1)+90<0
            WD_backup3(i,1) = WD_backup3(i,1)*(-1)+450;
        else
            WD_backup3(i,1) = WD_backup3(i,1)*(-1)+90;
        end
    end
    for i = 1:sizeb
        if data_backup3(i,1) > temp_max
            data_backup3(i,1) =temp_max;
        elseif data_backup3(i,1) < temp_min
            data_backup3(i,1) =temp_min;
        end
    end

    WD_backup4 = WD4;
    data_backup4 = data4;
    a = size(WD_backup4);
    sizeb = a(1,1);
    clear a
    for i = 1:sizeb
        if WD_backup4(i,1)*(-1)+90<0
            WD_backup4(i,1) = WD_backup4(i,1)*(-1)+450;
        else
            WD_backup4(i,1) = WD_backup4(i,1)*(-1)+90;
        end
    end
    for i = 1:sizeb
        if data_backup4(i,1) > temp_max
            data_backup4(i,1) =temp_max;
        elseif data_backup4(i,1) < temp_min
            data_backup4(i,1) =temp_min;
        end
    end

    figure((jw-1)*4+1)
%     subplot(1,4,1)
    wind_rose(WD_backup(:,1),data_backup(:,1),'ci',[2,4,6,8,10],'di',[0,2,4,6,8])

    figure((jw-1)*4+2)
% subplot(1,4,2)
    wind_rose(WD_backup2(:,1),data_backup2(:,1),'ci',[2,4,6,8,10],'di',[0,2,4,6,8])

    figure((jw-1)*4+3)
% subplot(1,4,3)
    wind_rose(WD_backup3(:,1),data_backup3(:,1),'ci',[2,4,6,8,10],'di',[0,2,4,6,8])

    figure((jw-1)*4+4)
% subplot(1,4,4)
    wind_rose(WD_backup4(:,1),data_backup4(:,1),'ci',[2,4,6,8,10],'di',[0,2,4,6,8])
jw
% R = input('=')    
end
clear jw


%%
clear sizeb
clear ans temp_max temp_min i


