function sys_info = PrepareWorkspace(params)

%% Set up timer
current_time = datetime('now', 'Format','yyyy-MM-dd_HH-mm');
timer.total.start = tic;

%% Define paths
% Current path
dir_path.current = cd;

% Define path to functions
dir_path.functions = {'functions', 'functions_auxiliary'};

% Set up path to functions
if isfield(dir_path, 'functions')
    for idx_dir = 1:length(dir_path.functions)
        dir_path.add2path = ...
            sprintf('%s/%s', ...
            dir_path.current, dir_path.functions{idx_dir});
        if exist(dir_path.add2path, 'dir')
            addpath(dir_path.add2path)
        end
    end
end

% Set up path to results
dir_path.results.main = sprintf('%s/%s', dir_path.current, 'results');
if ~exist(dir_path.results.main, 'dir')
    mkdir(dir_path.results.main)
end

% Subdirectory
dir_path.results.subdir = ...
    sprintf('%s/%s', dir_path.results.main, params.calc.model_name);
if isfield(params.calc, 'long_calc') && (params.calc.long_calc == true)
    dir_path.results.subdir = ...
        sprintf('%s_%s', dir_path.results.subdir, 'long_calc');
else
    dir_path.results.subdir = ...
        sprintf('%s_%s', dir_path.results.subdir, char(current_time));
end
make_dir(dir_path.results.subdir);

% Subdirectory and file for parameters
dir_path.params = dir_path.results.subdir;
make_dir(dir_path.params);
dir_path.params_file = ...
    sprintf('%s/parameters.txt', dir_path.results.subdir);

% Subdirectory for meshes
dir_path.mesh =  sprintf('%s/meshes', dir_path.results.subdir);
make_dir(dir_path.mesh)

% Subdirectory for geometry
dir_path.geom =  sprintf('%s/geometry', dir_path.results.subdir);
make_dir(dir_path.geom)

% Subdirectory for spectra
dir_path.spectra =  sprintf('%s/spectra', dir_path.results.subdir);
make_dir(dir_path.spectra)

% Subdirectory for band diagram data and plots
dir_path.spectra_data =  sprintf('%s', dir_path.results.subdir);
make_dir(dir_path.spectra_data)
% dir_path.spectra_plot = ...
%     sprintf('%s/spectra_plot', dir_path.results.subdir);
% make_dir(dir_path.spectra_plot)

%% Load style for figures
plot_style_ITMO


%% Index to show text in command window
cwtext.idx_row = 0;

sys_info.cwtext = cwtext;
sys_info.current_time = current_time;
sys_info.dir_path = dir_path;
sys_info.timer = timer;