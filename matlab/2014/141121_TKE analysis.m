%% History
% 141121 Mr.Je-Woo Hong for TKE budget analysis (with Junhong Lee and
% Prof.Jinkyu Hong)


%% main run % using "temp" from "extract" code
u_bar = mean(temp(:,1));
v_bar = mean(temp(:,2));
w_bar = mean(temp(:,3));
Ts_bar = mean(temp(:,4));
TKE = 0;
uw = 0;
vw = 0;
wTs = 0;
for i = 1:36000
    TKE = TKE + 1/2*((temp(i,1)-u_bar)^2 + (temp(i,2)-v_bar)^2 + (temp(i,3)-w_bar)^2 );
    uw = uw + (temp(i,1)-u_bar)*(temp(i,3)-w_bar);
    vw = vw + (temp(i,2)-v_bar)*(temp(i,3)-w_bar);
    wTs = wTs + (temp(i,3)-w_bar)*(temp(i,4)-Ts_bar);
end
clear i 
TKE = TKE/36000;
uw = uw/36000;
vw = vw/36000;
wTs = wTs/36000;


%%
result_TKE(1,1) = u_bar;
clear u_bar
result_TKE(1,2) = v_bar;
clear v_bar
result_TKE(1,3) = w_bar;
clear w_bar
result_TKE(1,4) = Ts_bar;
clear Ts_bar
result_TKE(1,5) = TKE;
clear TKE
result_TKE(1,6) = uw;
clear uw
result_TKE(1,7) = vw;
clear vw
result_TKE(1,8) = wTs;
clear wTs




