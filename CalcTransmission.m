function [cost_val, T0, T1 ,T2, Tm1, Tm2, R0, R1, R2, Rm1, Rm2, Ttot, Rtot, Stot, sys_info] = ...
    CalcTransmission(model, params, sys_info)

timer = sys_info.timer;
cwtext = sys_info.cwtext;
current_time = sys_info.current_time;
dir_path = sys_info.dir_path;

% Assign value of k to COMSOL variable
model.param.set('freq0', [num2str(params.wave.freq0), ' [Hz]']);

% Run calculations in COMSOL
model.study('std1').run;

% Calculate diffraction coefficients
T0 = mphglobal(model, {'T_0'});
T1 = mphglobal(model, {'T_1'});
T2 = mphglobal(model, {'T_2'});
Tm1 = mphglobal(model, {'T_m1'});
Tm2 = mphglobal(model, {'T_m2'});
R0 = mphglobal(model, {'R_0'});
R1 = mphglobal(model, {'R_1'});
R2 = mphglobal(model, {'R_2'});
Rm1 = mphglobal(model, {'R_m1'});
Rm2 = mphglobal(model, {'R_m2'});
Ttot = mphglobal(model, {'T_tot'});
Rtot = mphglobal(model, {'R_tot'});
Stot = mphglobal(model, {'scat_tot'});


if Ttot > 1.1 || Rtot > 1.1 || Stot > 1.1
    cost_val = inf;
else
%     cost_val = min([abs(1 - Tm1), abs(1 - T1)]);
%     cost_val = min([abs(1 - Tm1), abs(1 - T1)]);
cost_val = abs(1 - Tm1);
end

%Save the results
dir_path.data_Tavg = ...
    sprintf('%s/Tavg.txt', dir_path.results.subdir);
dir_path.file2write_Tavg = fopen(dir_path.data_Tavg, 'a+');
formatspec = ['%i', repmat(', %.10f', 1, params.elements.n_vertices*2 + 3), '\n'];

vertices_coords = ...
    reshape([cell2mat(params.elements.vx); cell2mat(params.elements.vy)], 1,[]);

array2save = [sys_info.calc_tag, cost_val, T1, Tm1, vertices_coords];
fprintf(dir_path.file2write_Tavg, formatspec, array2save');

fclose('all');

sys_info.timer = timer;
sys_info.cwtext = cwtext;
sys_info.current_time = current_time;
sys_info.dir_path = dir_path;