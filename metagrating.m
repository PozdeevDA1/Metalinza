function z = metagrating(x, params, sys_info)

cwtext = sys_info.cwtext;
current_time = sys_info.current_time;
dir_path = sys_info.dir_path;
timer = sys_info.timer;

% %% Show info in the command line:
% cwtext.idx_row = cwtext.idx_row + 1;
% cwtext.row{cwtext.idx_row} =...
%     sprintf('Evaluated defects: 0 from %i \n', params.algorithm.n_evaluations);
% print_cwtext(cwtext)
% 
% for idx_defect = 1:params.algorithm.n_evaluations
%     n_resonators2remove = ...
%         randi([n_resonators2remove_min, n_resonators2remove_max]);
%     resonators2remove = randi([1, n_resonators], 1, n_resonators2remove);

params.elements.remove = find(x == 1);



% Vertices:
params.elements.vx = num2cell(x(1:2:end));
params.elements.vy = num2cell(x(2:2:end));

params.unit_cell.ay = max(cell2mat(params.elements.vy)) - ...
    min(cell2mat(params.elements.vy));

% Array of elements: geometric size
params.elements.width = max(cell2mat(params.elements.vx)) - ...
    min(cell2mat(params.elements.vx));

params.elements.height = abs(max(cell2mat(params.elements.vy)));

% Host

params.host.size_y = 3*params.wave.ff_zone + params.elements.height;


% Integration lines
params.host.R_line_y = params.elements.height + params.wave.ff_zone;
params.host.T_line_y = (min(cell2mat(params.elements.vy)))-params.wave.ff_zone;

% model = SetupModel(params, sys_info);
% 
% [cost_val, ~, ~ ,~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
%     CalcTransmission(model, params, sys_info);
try
model = SetupModel(params, sys_info);

[cost_val, ~, ~ ,~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    CalcTransmission(model, params, sys_info);

catch
    cost_val = inf;
    %Save the results
    dir_path.data_Tavg = ...
    sprintf('%s/Tavg.txt', dir_path.results.subdir);
    dir_path.file2write_Tavg = fopen(dir_path.data_Tavg, 'a+');
    formatspec = ['%i', repmat(', %.10f', 1, params.elements.n_vertices*2 + 3), '\n'];

    vertices_coords = ...
        reshape([cell2mat(params.elements.vx); cell2mat(params.elements.vy)], 1,[]);

    T1 = 0;
    Tm1 = 0;
    array2save = [sys_info.calc_tag, cost_val, T1, Tm1, vertices_coords];
    fprintf(dir_path.file2write_Tavg, formatspec, array2save');

end

%     cwtext.row{cwtext.idx_row} =...
%         sprintf('Evaluated defects: %i from %i \n', ...
%         idx_defect, params.algorithm.n_evaluations);
% end

z = cost_val;
