%% calculation
% data (1: elevation, 2: wind-direction, 3: albedo)

result = zeros(length(data),2); % 1: distance, 2: wind-direction
result(:,2) = data(:,2);
for i = 1:length(data)
    result(i,1) = 27.7/tan(data(i,1)/180*pi());
end
clear i

result_xy = zeros(length(data),2); % 1: x-east, 2: y-north
for i = 1:length(data)
    if result(i,2)<90
        result_xy(i,1) = result(i,1)*cos((90-result(i,2))/180*pi());
        result_xy(i,2) = result(i,1)*sin((90-result(i,2))/180*pi());
    elseif result(i,2)<180
        result_xy(i,1) = result(i,1)*cos((result(i,2)-90)/180*pi());        
        result_xy(i,2) = -1*result(i,1)*sin((result(i,2)-90)/180*pi());
    elseif result(i,2)<270
        result_xy(i,1) = -1*result(i,1)*cos((270-result(i,2))/180*pi());
        result_xy(i,2) = -1*result(i,1)*sin((270-result(i,2))/180*pi());        
    else %result(i,2)>270
        result_xy(i,1) = -1*result(i,1)*cos((result(i,2)-270)/180*pi());
        result_xy(i,2) = result(i,1)*sin((result(i,2)-270)/180*pi());                
    end
end
clear i


%% figure
pcolor(dem);figure(gcf);
hold on
plot(200,207,'x','LineWidth',2,'MarkerSize',10);
result_xy(:,3) = 200; % radius
scatter(200+result_xy(:,1),207-result_xy(:,2),result_xy(:,3),data(:,3),'filled');