%clear all;clc;


% Define your FRED API key (replace 'your_api_key' with your actual key)
apiKey = '1729e2c2e3908c169a742b6e953a69e7';

% Define the FRED API URL and series ID

baseUrl = 'https://api.stlouisfed.org/fred/series/observations';

seriesID = {'FEDFUNDS','GDPC1','GDPPOT','PCEPI'}; % Note: Make sure to use the correct series ID
seriesF  = {'M','Q','Q','M'}; % Note: FRED frequency of seriesID before aggregation

% first obs =
first_obs_q = '1970-01-01'; % 1970:Q1
first_obs_m = '1970-01-01'; % 1970:M1
% last obs
last_obs_q  = '2008-10-01'; % 2008:Q4
last_obs_m  = '2008-12-01'; % 2008:M12


% Initialize a structure to store data
rawdata  = struct();
datam     = [];
dataq     = [];
dates    = struct(); % Structure to store first and last dates for each series

% Loop through each series ID
for i = 1:length(seriesID)
    % Construct the full API request URL with the frequency parameter
    url = sprintf('%s?series_id=%s&api_key=%s&file_type=json&value_format=full', ...
                  baseUrl, seriesID{i}, apiKey);

    % Fetch the data
    response = webread(url);

    % Store the response in a field named after the series ID
    rawdata.(seriesID{i})  = str2double(cellstr(char(response.observations.value)));% response;
    dates.(seriesID{i})    = char(response.observations.date);

    % % Extract the first and last observation dates
    
    switch seriesF{i}
        case 'Q'
            first_obs_index=find(all(dates.(seriesID{i})==first_obs_q,2)==1);
            last_obs_index =find(all(dates.(seriesID{i})==last_obs_q,2)==1);

 
            dataq    = [dataq,rawdata.(seriesID{i})(first_obs_index:last_obs_index,1)];

        case 'M'
            first_obs_index=find(all(dates.(seriesID{i})==first_obs_m,2)==1);
            last_obs_index =find(all(dates.(seriesID{i})==last_obs_m,2)==1);

            %datam     = [datam,rawdata.(seriesID{i})(first_obs_index:last_obs_index,1)];
            datesm       = dates.(seriesID{i})(first_obs_index:last_obs_index,:);
            datetime_vec = datetime(datesm,'InputFormat','yyyy-MM-dd'); 

            TTm = timetable(datetime_vec,rawdata.(seriesID{i})(first_obs_index:last_obs_index,1));
            TTq = convert2quarterly(TTm,'Aggregation','mean');

            dataq = [dataq,TTq.Var1];
            datesq.(seriesID{i}) = datestr(TTq.datetime_vec,'yyyy-mm-dd'); 
 
    end
    % observations = response.observations;
    % if ~isempty(observations)
    %     firstDate = observations(1).date;
    % 
    % 
    %     lastDate  = observations(end).date;
    % 
    % 
    %     dates.(seriesID{i}) = struct('FirstDate', firstDate, 'LastDate', lastDate);
    % else
    %     dates.(seriesID{i}) = struct('FirstDate', 'N/A', 'LastDate', 'N/A');
    % end
end

% aggregate monthly data to quarterly

