function Save_Data_mini2p(data, activityMap, valid_PCs)

data.activityMap = activityMap;
data.PCs = valid_PCs;

[y,m,d] = ymd(datetime("today"));
y = num2str(y);
m = num2str(m);
d = num2str(d);
if length(d) == 1 % if it's a one-digit day
    d = strcat('0',d);
end
if length(m) == 1
    m = strcat('0',m);
end
save(strcat(y(3:end),m,d,'_processed_data_MoserCriteria.mat'), 'data'); % saves to file with current date (YYMMDD)
disp('Done! Processed data saved to current directory.')