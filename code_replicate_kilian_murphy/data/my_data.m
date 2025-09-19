clear all;clc;



% Define your FRED API key (replace 'your_api_key' with your actual key)
apiKey = '1729e2c2e3908c169a742b6e953a69e7';

% Define the FRED API URL and series ID

baseUrl = 'https://api.stlouisfed.org/fred/series/observations';

seriesID = {'FEDFUNDS','GDPC1','GDPPOT','GDPDEF'}; % Note: Make sure to use the correct series ID


% Initialize a structure to store data
data = struct();
dates = struct(); % Structure to store first and last dates for each series

% Loop through each series ID
for i = 1:length(seriesID)
    % Construct the full API request URL with the frequency parameter
    url = sprintf('%s?series_id=%s&api_key=%s&file_type=json&frequency=q', ...
                  baseUrl, seriesID{i}, apiKey);

    % Fetch the data
    response = webread(url);

    % Store the response in a field named after the series ID
    data.(seriesID{i}) = response;

    % Extract the first and last observation dates
    observations = response.observations;
    if ~isempty(observations)
        firstDate = observations(1).date;
      
       
        lastDate  = observations(end).date;


        dates.(seriesID{i}) = struct('FirstDate', firstDate, 'LastDate', lastDate);
    else
        dates.(seriesID{i}) = struct('FirstDate', 'N/A', 'LastDate', 'N/A');
    end
end