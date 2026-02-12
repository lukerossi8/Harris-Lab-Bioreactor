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
    filepath = "../Bioreactor Simulation Flows/Bioreactor_data_7deg_20rpm_lv6_onecycle/" + filename;

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

    all_data(:,:,i-2) = file_2d_array;

    % Extracting the time from the file name
    to_remove = ["Data_all_64_", "_0.dat"];
    time_str = erase(filename, to_remove);
    time = str2double(time_str);
    times(i-2) = time;
end

snapshot1 = all_data(:,:,67);
snapshot2 = all_data(:,:,4);

% Dimensionalizing
% Reference frames
L_ref = 0.25; % m
L_y = 1/3.5; % m
H_ref = L_ref*L_y; % m2 I guess

theta_max = 7*pi/180; % in radians
rpm = 20; % in rpm, given in flow data
T_per = 60/rpm; % period in seconds

V_ref = L_ref/4*(H_ref + 0.5*L_ref*tan(theta_max));
U_ref = V_ref/(H_ref*0.5)/T_per; % characteristic velocity scale in m/s I guess


%  H_bio  = L_bio*Ly;
%   V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
%   U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity scale
% 1:40
% L_bio = 0.25m
% 1:40
% Ly = 0.286 = *(1/3.5)

save("time_dependent_flow_data_out.mat", "all_data", "times");

