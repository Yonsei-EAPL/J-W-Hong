%% information for the position

% for Eunpyeong Newtown
i1 = 25;
j1 = 22;

% for Seoul Forest
i2 = 20;
j2 = 24;


%% call the data

rainnc = nc_varget('sum_2','RAINNC');
rainc = nc_varget('sum_2','RAINC');
tsk = nc_varget('sum_2','TSK');
t2 = nc_varget('sum_2','T2');
u10 = nc_varget('sum_2','U10');
v10 = nc_varget('sum_2','V10');
swdown = nc_varget('sum_2','SWDOWN');
glw = nc_varget('sum_2','GLW');
hfx = nc_varget('sum_2','HFX');
lh = nc_varget('sum_2','LH');


%% variable for Eunpyeong Newtown and Seoul Forest

% for Eunpyeong Newtown
E_rain = zeros(2947,1); % rain total (mm)
E_ts = zeros(2947,1); % surface temp (K)
E_t2 = zeros(2947,1); % 2m temp (K)
E_u10 = zeros(2947,1); % 10m u (ms-1)
E_v10 = zeros(2947,1); % 10m v (ms-1)
E_swdn = zeros(2947,1); % sw down (Wm-2)
E_lwdn = zeros(2947,1); % lw down (Wm-2)
E_h = zeros(2947,1); % sensible heat flux (Wm-2)
E_le = zeros(2947,1); % latent heat flux (Wm-2)

% for Seoul Forest
S_rain = zeros(2947,1);
S_ts = zeros(2947,1); 
S_t2 = zeros(2947,1);
S_u10 = zeros(2947,1);
S_v10 = zeros(2947,1);
S_swdn = zeros(2947,1);
S_lwdn = zeros(2947,1);
S_h = zeros(2947,1);
S_le = zeros(2947,1);


%% extract the data

for i = 1:2947
    E_rain(i,1) = rainc(i,i1,j1) + rainnc(i,i1,j1);
    E_ts(i,1) = tsk(i,i1,j1);
    E_t2(i,1) = t2(i,i1,j1);
    E_u10(i,1) = u10(i,i1,j1);
    E_v10(i,1) = v10(i,i1,j1);
    E_swdn(i,1) = swdown(i,i1,j1);
    E_lwdn(i,1) = glw(i,i1,j1);
    E_h(i,1) = hfx(i,i1,j1);
    E_le(i,1) = lh(i,i1,j1);
    
    S_rain(i,1) = rainc(i,i2,j2) + rainnc(i,i2,j2);
    S_ts(i,1) = tsk(i,i2,j2);
    S_t2(i,1) = t2(i,i2,j2);
    S_u10(i,1) = u10(i,i2,j2);
    S_v10(i,1) = v10(i,i2,j2);
    S_swdn(i,1) = swdown(i,i2,j2);
    S_lwdn(i,1) = glw(i,i2,j2);
    S_h(i,1) = hfx(i,i2,j2);
    S_le(i,1) = lh(i,i2,j2);
end

clear i 


