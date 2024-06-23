%% Preparation 

close all;
clear;
clc;

%% Prepare workspace (paths, timer, etc.)
params.calc.model_name = 'polygons';
params.calc.algorithm = 'GA';
params.calc.hard_boundaries = true;

sys_info = PrepareWorkspace(params);
timer = sys_info.timer;
cwtext = sys_info.cwtext;
current_time = sys_info.current_time;
dir_path = sys_info.dir_path;

%% Load parameters
params = LoadParameters(params);

%% Problem definition
problem.CostFun = @(x) metagrating(x, params, sys_info);

problem.n_var = params.elements.n_vertices;

%% Show information in command window
cwtext.idx_row = cwtext.idx_row + 1;
cwtext.row{cwtext.idx_row} =...
    sprintf([...
    'Transmission of a resonators wall with defects.\n', ...
    'Calculations started at %s.\n\n'], ...
    char(current_time));
print_cwtext(cwtext)


%% Save parameters
params_title = sprintf(['Model: %s \n', 'Launched: %s\n'], ...
    params.calc.model_name, char(current_time));
file2write_params = fopen(dir_path.params_file, 'a+');
param_file_content = dir(dir_path.params_file);
if param_file_content.bytes == 0
    save_params(params, dir_path.params_file, params_title);
end


%% Algorithm implementation

problem.n_var = params.elements.n_vertices*2;

vx_min = -params.unit_cell.period_x/2 + params.unit_cell.min_element;
vx_max = params.unit_cell.period_x/2 - params.unit_cell.min_element;
vy_min = 0;
vy_max = 2*params.unit_cell.period_x;

problem.VarMin = ...
    reshape([repmat(vx_min, 1, params.elements.n_vertices);...
    repmat(vy_min, 1, params.elements.n_vertices)], 1,[]);

problem.VarMax = reshape([repmat(vx_max, 1, params.elements.n_vertices);...
    repmat(vy_max, 1, params.elements.n_vertices)], 1,[]);


title_str = '';
for idx_v = 1:params.elements.n_vertices
    title_str = [title_str, sprintf(', x_%i, y_%i', idx_v, idx_v)];
end
dir_path.data_Tavg = ...
    sprintf('%s/Tavg.txt', dir_path.results.subdir);
dir_path.file2write_Tavg = fopen(dir_path.data_Tavg, 'a+');
fprintf(dir_path.file2write_Tavg, ...
    sprintf('N_iteration, Cost, T1, Tm1%s \n\n', title_str));

sys_info.dir_path = dir_path;
% Run GA
output = RealGeneticAlgorithm(problem, params, sys_info);


%% Finalization
% Calculation time
timer.total.end = toc(timer.total.start);
timer.total.end = format_time(timer.total.end);

% Show information in command window
cwtext.idx_row = cwtext.idx_row + 1;
cwtext.row{cwtext.idx_row} = sprintf(['\n', ...
    'Calculations are completed.\n', ...
    'Elapsed time is %s\n'], timer.total.end);
print_cwtext(cwtext)

%% Post-processig
figure;
plot(1:params.algorithm.max_iter, output.best_cost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;