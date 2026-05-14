clc; clearvars;

%% Constructing a directory of all input data
% Each file contains the flow data at a given timestep
files = struct2table(dir("**")); % Creating a directory of all files in the current folder and all subfolders
data_folder = "Bioreactor_data_7deg_20rpm_lv6_onecycle"; % must be a direct subfolder of current folder
in_data_folder = contains(files.folder, data_folder); % Identifying the files in the desired subfolder
data_files = files(in_data_folder, :); % Filtering for only the files in the desired subfolder, i.e. the data files

all_data = []; % Initializing array to hold all the data
times = zeros(height(data_files)-2, 1); % Initializing array to hold each timestep associated with a data file

for i=3:height(data_files) % Ignoring first two files, '.' and '..'
    
    filename = data_files{i, "name"};
    filepath = "./Bioreactor_data_7deg_20rpm_lv6_onecycle/" + filename;

    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 13);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = " ";

    % Specify column names and types
    opts.VariableNames = ["x", "y", "ux", "uy", "vol_frac", "ufx", "ufy", "cs", "fsx", "fsy", "fmx", "fmy", "vorticity"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";

    % Import the data
    file_table_raw = readtable(filepath, opts);
    file_table_sort1 = sortrows(file_table_raw,1);
    file_table_sort2 = sortrows(file_table_sort1,2,'descend');
    file_2d_array = table2array(file_table_sort2); % This is the sorted 2d array of the data from the particular file
    file_2d_array_trim = file_2d_array(1409:2688, :); % Manually extracting the non-empty rows, as many only had blank values


    all_data(:,:,i-2) = file_2d_array_trim;

    % Extracting the time from the file name
    to_remove = ["Data_all_64_", "_0.dat"];
    time_str = erase(filename, to_remove);
    time = str2double(time_str);
    times(i-2) = time;
end

% Normalizing times to start at 0
times = times - times(1);

% Dimensionalizing
L_ref = 0.25; % characteristic length scale, m (bioreactor length)
T_per = 3; % s (period)
non_dim_per = times(end) + mean(diff(times)); % non-dimensionalized period
T_ref = T_per/non_dim_per; % characteristic time scale, s
V_ref = L_ref/T_ref; % characteristic velocity scale, m/s

times = times.*T_ref;

x = all_data(:,1,:);
y = all_data(:,2,:);
vx = all_data(:,3,:);
vy = all_data(:,4,:);

x_dim = x.*L_ref; % m
y_dim = y.*L_ref; % m
vx_dim = vx.*V_ref; % m/s
vy_dim = vy.*V_ref; % m/s

all_data(:,1,:) = x_dim;
all_data(:,2,:) = y_dim;
all_data(:,3,:) = vx_dim;
all_data(:,4,:) = vy_dim;

save("time_dependent_flow_data_out.mat", "all_data", "times");

