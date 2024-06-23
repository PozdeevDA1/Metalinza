function model = SetupModel(params, sys_info)

dir_path = sys_info.dir_path;

%% Create a model
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
comp1 = model.component.create('comp1');

%% Default parameters
% model.param.set('freq0', [num2str(params.wave.freq0_start), ' [Hz]']);

%% Create the geometry
geom0 = comp1.geom.create('geom0',2);

%% Host
host_coord_x = -params.unit_cell.period_x/2;%-(params.host.size_x - params.elements.width)/2;
host_coord_y = -params.host.size_y/2 + params.elements.height/2;
host = geom0.feature.create('host','Rectangle');
host.set('size', [params.host.size_x, params.host.size_y]);
host.set('pos', [host_coord_x, host_coord_y]);
host.set('base', 'corner');
host.set('createselection', 'on');

%% PML
PML_bottom = geom0.feature.create('PML_bottom','Rectangle');
PML_bottom.set('size', ...
    [params.host.size_x, params.pml.width]);
PML_bottom_coord_x = host_coord_x;
PML_bottom_coord_y = host_coord_y - params.pml.width;
PML_bottom.set('pos', [PML_bottom_coord_x, PML_bottom_coord_y]);
PML_bottom.set('base', 'corner');
PML_bottom.set('createselection', 'on');

PML_top = geom0.feature.create('PML_top','Rectangle');
PML_top.set('size', ...
    [params.host.size_x, params.pml.width]);
PML_top_coord_x = host_coord_x;
PML_top_coord_y = host_coord_y + params.host.size_y;
PML_top.set('pos', [PML_top_coord_x, PML_top_coord_y]);
PML_top.set('base', 'corner');
PML_top.set('createselection', 'on');

%% Sampling lines
% Reflection
ls1_coord_1 = {num2str(host_coord_x) num2str(params.host.R_line_y)};
ls1_coord_2 = {num2str(host_coord_x + params.host.size_x) num2str(params.host.R_line_y)};
ls1 = geom0.create('ls1', 'LineSegment');
ls1.set('specify1', 'coord');
ls1.set('coord1', ls1_coord_1);
ls1.set('specify2', 'coord');
ls1.set('coord2', ls1_coord_2);

% Transmission
ls2_coord_1 = {num2str(host_coord_x) num2str(params.host.T_line_y)};
ls2_coord_2 = {num2str(host_coord_x + params.host.size_x) num2str(params.host.T_line_y)};
ls2 = geom0.create('ls2', 'LineSegment');
ls2.set('specify1', 'coord');
ls2.set('coord1', ls2_coord_1);
ls2.set('specify2', 'coord');
ls2.set('coord2', ls2_coord_2);

%% Array of the elements
label = 0;
label = label + 1;
polygon_tag = sprintf('element_%i', label);
poly1 = geom0.feature.create(polygon_tag, 'Polygon');
poly1.set('source', 'table');

for idx_v = 1:params.elements.n_vertices
    poly1.setIndex('table', params.elements.vx{idx_v}, idx_v - 1, 0);
    poly1.setIndex('table', params.elements.vy{idx_v}, idx_v - 1, 1);
end
poly1.setIndex('table', params.elements.vx{1}, 0, 0);
poly1.setIndex('table', params.elements.vy{1}, 0, 1);
poly1.set('createselection', 'on');


% Save geometry
dir_path.fig_geom = ...
    sprintf('%s/%i_geom.png', dir_path.geom, sys_info.calc_tag);
fig_geom = figure('visible','off');
mphgeom(model,'geom0');
xlabel('x')
ylabel('y')
ax = gca;
% ax.YLim = [-0.5*params.unit_cell.ay 1.5*params.unit_cell.ay];

try
    exportgraphics(gcf, dir_path.fig_geom, 'Resolution',300);
catch
    saveas(gcf, dir_path.fig_geom);
end
close(fig_geom)
clear fig_geom

%% Set up domain labels and selections
% PML
PML_position = {'top', 'bottom'};

PML_domains = [];
for idx_position = 1:length(PML_position)
    PML_domains_ent = ...
        mphgetselection(model.selection(...
        sprintf('geom0_PML_%s_dom', PML_position{idx_position})));
    PML_domains = [PML_domains, PML_domains_ent.entities];

    PML_boundaries_ent = ...
        mphgetselection(model.selection(...
        sprintf('geom0_PML_%s_bnd', PML_position{idx_position})));
    PML_boundaries.(PML_position{idx_position}) = ...
        PML_boundaries_ent.entities;
end
comp1.selection.create('sel_PML', 'Explicit');
comp1.selection('sel_PML').set(PML_domains);
comp1.selection('sel_PML').label('PML');

PML_bottom_left = ...
    mphselectbox(model, 'geom0',...
    [PML_bottom_coord_x - 0.5*params.host.size_x, ...
    PML_bottom_coord_y - 0.1*params.pml.width; ...
    PML_bottom_coord_x + 0.5*params.host.size_x, ...
    PML_bottom_coord_y + 1.1*params.pml.width]', 'boundary');

PML_bottom_right = ...
    mphselectbox(model, 'geom0',...
    [PML_bottom_coord_x + 0.5*params.host.size_x, ...
    PML_bottom_coord_y - 0.1*params.pml.width; ...
    PML_bottom_coord_x + 1.1*params.host.size_x, ...
    PML_bottom_coord_y + 1.1*params.pml.width]', ...
    'boundary');

PML_top_left = ...
    mphselectbox(model, 'geom0',...
    [PML_top_coord_x - 0.5*params.host.size_x, ...
    PML_top_coord_y - 0.1*params.pml.width; ...
    PML_top_coord_x + 0.5*params.host.size_x, ...
    PML_top_coord_y + 1.1*params.pml.width]', 'boundary');

PML_top_right = ...
    mphselectbox(model, 'geom0',...
    [PML_top_coord_x + 0.5*params.host.size_x, ...
    PML_top_coord_y - 0.1*params.pml.width; ...
    PML_top_coord_x + 1.1*params.host.size_x, ...
    PML_top_coord_y + 1.1*params.pml.width]', ...
    'boundary');
PML_boundaries_mesh = ...
    [PML_bottom_left, PML_bottom_right, PML_top_left, PML_top_right];


% % Elements
elements_domains = [];
params.elements.n_elements = 1;
elements_domains_ent = ...
    mphgetselection(model.selection(...
    sprintf('geom0_element_%i_dom', 1)));
elements_domains = ...
    [elements_domains, elements_domains_ent.entities];

comp1.selection.create('sel_elements', 'Explicit');
comp1.selection('sel_elements').set(elements_domains);
comp1.selection('sel_elements').label('elements');

% Host
host_domain_ent = mphgetselection(model.selection('geom0_host_dom'));
host_domain_bug = host_domain_ent.entities;
host_domain = setdiff(host_domain_bug, elements_domains);
comp1.selection.create('sel_host', 'Explicit');
comp1.selection('sel_host').set([host_domain, PML_domains]);
comp1.selection('sel_host').label('host');

% model.selection.tags
% mphgetselection(model.selection('geom0_host_dom'))

% Host boundaries
coordbox_left = [host_coord_x - 0.1*host_coord_x, host_coord_y - 1.1*params.pml.width; ...
    host_coord_x + 0.9*params.unit_cell.min_element, host_coord_y + params.host.size_y + 1.5*params.pml.width];
host_left = ...
    mphselectbox(model, 'geom0', coordbox_left', 'boundary');

coordbox_right = [host_coord_x + params.host.size_x - 0.9*params.unit_cell.min_element, host_coord_y - 1.1*params.pml.width; ...
    host_coord_x + 1.1*params.host.size_x, host_coord_y + params.host.size_y + 1.5*params.pml.width];
host_right = ...
    mphselectbox(model, 'geom0', coordbox_right', 'boundary');


%% Set up physics
acpr_phys = comp1.physics.create('acpr', 'PressureAcoustics', 'geom1');
acpr_phys.selection.set([host_domain, PML_domains]);
acpr_phys.prop('ReferencePressure').set('ReferenceType', 'ReferencePressureWater');
acpr_phys.prop('cref').set('cref', params.wave.c_water);

acpr_phys.feature.create('bpf1', 'BackgroundPressureField', 2);
acpr_phys.feature('bpf1').selection.set(host_domain);
acpr_phys.feature('bpf1').set('PressureFieldType', 'UserDefined');
acpr_phys.feature('bpf1').set('p', 'p_inc');

floquet = cell(1, length(host_left));
for idx_bound = 1:length(host_left)
    floquet{idx_bound} = acpr_phys.feature.create(sprintf('pc%i', idx_bound), 'PeriodicCondition', 1);
    floquet{idx_bound}.set('PeriodicType', 'Floquet');
    floquet{idx_bound}.set('kFloquet', {'kx'; 'ky'; '0'});
    floquet{idx_bound}.selection.set([host_left(idx_bound), host_right(idx_bound)]);
end

if ~params.calc.hard_boundaries
    solid_phys = comp1.physics.create('solid', 'SolidMechanics', 'geom0');
    solid_phys.selection.set(elements_domains);

    comp1.multiphysics.create('asb1', 'AcousticStructureBoundary', 1);
    comp1.multiphysics('asb1').selection.all;
end

%% Materials
mat_water = model.component('comp1').material.create('mat_water', 'Common');
load_water
mat_water.selection.set([host_domain, PML_domains]);
% 
% mat_aluminium = model.component('comp1').material.create('mat_aluminium', 'Common');
% load_aluminium
% mat_aluminium.selection.set(elements_domains);

mat_plastic = model.component('comp1').material.create('mat_plastic', 'Common');
load_plastic
mat_plastic.selection.set(elements_domains);


%% Set up PML
comp1.coordSystem.create('pml1', 'PML');
comp1.coordSystem('pml1').selection.set(PML_domains);
comp1.coordSystem('pml1').set('PMLfactor', num2str(params.pml.scaling));
comp1.coordSystem('pml1').set('PMLgamma', num2str(params.pml.curvature));


%% Mesh
comp1.mesh.create('mesh1');

% Size params
comp1.mesh('mesh1').feature('size').set('hmax', params.mesh.max_element);
comp1.mesh('mesh1').feature('size').set('hmin', params.mesh.min_element);
comp1.mesh('mesh1').feature('size').set('hgrad', params.mesh.growth_rate);
comp1.mesh('mesh1').feature('size').set('hcurve', params.mesh.curvature);
comp1.mesh('mesh1').feature('size').set('hnarrow', params.mesh.narrow_resolution);

% PML
comp1.mesh('mesh1').feature('size').set('custom', true);
comp1.mesh('mesh1').create('map1', 'Map');
comp1.mesh('mesh1').feature('map1').selection.geom('geom0', 2);
comp1.mesh('mesh1').feature('map1').selection.named('sel_PML');
comp1.mesh('mesh1').feature('map1').create('dis1', 'Distribution');
comp1.mesh('mesh1').feature('map1').feature('dis1').selection.set(PML_boundaries_mesh);
comp1.mesh('mesh1').feature('map1').feature('dis1').set(...
    'numelem', params.mesh.pml_distribution);

% Host
comp1.mesh('mesh1').create('ftri1', 'FreeTri');
comp1.mesh('mesh1').feature('ftri1').selection.remaining;
% comp1.mesh('mesh1').feature('ftri1').selection.geom('geom0', 2);
% mphsave(model,sprintf('%s/%i_model', dir_path.geom, sys_info.calc_tag),'component','on')
% comp1.mesh('mesh1').feature('ftri1').selection.set(host_domain);
% 
% % Elements
% comp1.mesh('mesh1').create('ftri2', 'FreeTri');
% comp1.mesh('mesh1').feature('ftri2').selection.geom('geom0', 2);
% % comp1.mesh('mesh1').feature('ftri2').selection.set(elements_domains);
% comp1.mesh('mesh1').feature('ftri2').selection.remaining;
% mphsave(model,sprintf('%s/%i_model', dir_path.geom, sys_info.calc_tag),'component','on');
comp1.mesh('mesh1').run;

% Save mesh
if params.mesh.save_mesh
    dir_path.fig_mesh = ...
        sprintf('%s/%i_mesh.png', dir_path.mesh, sys_info.calc_tag);
    fig_mesh = figure('visible','off');
    mphmesh(model);
    xlabel('x')
    ylabel('y')
    saveas(gcf, dir_path.fig_mesh);
    close(fig_mesh)
    clear fig_mesh
end

%% Set up study
model.study.create('std1');
model.study('std1').create('freq', 'Frequency');
model.study('std1').feature('freq').set('plist', params.wave.freq0);


%% Set up variables
% Integration operations
int_point_bnd = ...
    mphselectbox(model, 'geom0',...
    [PML_top_coord_x - 0.1*params.pml.width, ...
    PML_top_coord_y - 0.1*params.pml.width; ...
    PML_top_coord_x + 0.1*params.pml.width, ...
    PML_top_coord_y + 0.1*params.pml.width]', ...
    'point');
int_point = comp1.cpl.create('intop1', 'Integration');
int_point.set('axisym', true);
int_point.selection.geom('geom0', 0);
int_point.selection.set(int_point_bnd);

int_line1_bnd = ...
    mphselectbox(model, 'geom0',...
    [str2double(ls1_coord_1{1}) - params.unit_cell.min_element/2, ...
    str2double(ls1_coord_1{2}) - params.unit_cell.min_element; ...
    str2double(ls1_coord_2{1}) + params.unit_cell.min_element/2, ...
    str2double(ls1_coord_2{2}) + params.unit_cell.min_element]', 'boundary');
int_line1 = comp1.cpl.create('intop3', 'Integration');
int_line1.set('axisym', true);
int_line1.selection.geom('geom0', 1);
int_line1.selection.set(int_line1_bnd);

int_line2_bnd = ...
    mphselectbox(model, 'geom0',...
    [str2double(ls2_coord_1{1}) - params.unit_cell.min_element/2, ...
    str2double(ls2_coord_1{2}) - params.unit_cell.min_element; ...
    str2double(ls2_coord_2{1}) + params.unit_cell.min_element/2, ...
    str2double(ls2_coord_2{2}) + params.unit_cell.min_element]', 'boundary');
int_line2 = comp1.cpl.create('intop2', 'Integration');
int_line2.set('axisym', true);
int_line2.selection.geom('geom0', 1);
int_line2.selection.set(int_line2_bnd);

% Parameters
model.param.set('theta_inc', [num2str(params.wave.theta_inc), ' [deg]']);
model.param.set('p0', [num2str(params.wave.p0), ' [Pa]']);
model.param.set('period_x', [num2str(params.unit_cell.period_x), ' [m]']);
model.param.set('lam0', [num2str(params.wave.lam0), ' [m]']);

% Variables
var_inc_field = comp1.variable.create('var_inc_field');
var_inc_field.label('Variables Incident field');
var_inc_field.set('kx_e', 'sin(theta_inc)');
var_inc_field.set('ky_e', 'cos(theta_inc)');
var_inc_field.set('k0', 'intop1(acpr.k)');
var_inc_field.set('kx', 'k0*kx_e');
var_inc_field.set('ky', '-k0*ky_e');
var_inc_field.set('p_inc', 'p0*exp(-i*kx*x-i*ky*y)');
var_inc_field.set('A_inc', 'intop2(acpr.p_b*exp(i*kx*x+i*ky*y))/period_x');
var_inc_field.set('B_inc', 'intop3(acpr.p_b*exp(i*kx*x+i*ky*y))/period_x');

var_T0 = comp1.variable.create('var_T0');
var_T0.label('Variables T0');
var_T0.set('k_t_x_0', 'k0*sin(theta_inc)');
var_T0.set('k_t_y_0', '-sqrt(k0^2-k_t_x_0^2)');
var_T0.set('theta_t_0', 'theta_inc');
var_T0.set('A_0', 'intop2(acpr.p_t*exp(i*k_t_x_0*x+i*k_t_y_0*y))/period_x');
var_T0.set('T_0', '(abs(A_0))^2*cos(theta_t_0)/((abs(A_inc))^2*cos(theta_inc))');

var_T1 = comp1.variable.create('var_T1');
var_T1.label('Variables T1');
var_T1.set('k_t_x_1', 'k0*sin(theta_inc)-2*pi/period_x');
var_T1.set('k_t_y_1', '-sqrt(k0^2-k_t_x_1^2)');	
var_T1.set('theta_t_1', 'asin(sin(theta_inc)-lam0/period_x)');
var_T1.set('A_1', 'intop2(acpr.p_t*exp(i*k_t_x_1*x+i*k_t_y_1*y))/period_x');
var_T1.set('T_1', '(abs(A_1))^2*cos(theta_t_1)/((abs(A_inc))^2*cos(theta_inc))');

var_T2 = comp1.variable.create('var_T2');
var_T2.label('Variables T2');
var_T2.set('k_t_x_2', 'k0*sin(theta_inc)-4*pi/period_x');	
var_T2.set('k_t_y_2', '-sqrt(k0^2-k_t_x_2^2)');
var_T2.set('theta_t_2', 'asin(sin(theta_inc)-2*lam0/period_x)');
var_T2.set('A_2', 'intop2(acpr.p_t*exp(i*k_t_x_2*x+i*k_t_y_2*y))/period_x');	
var_T2.set('T_2', '(abs(A_2))^2*cos(theta_t_2)/((abs(A_inc))^2*cos(theta_inc))');

var_Tm1 = comp1.variable.create('var_Tm1');
var_Tm1.label('Variables Tm1');
var_Tm1.set('k_t_x_m1', 'k0*sin(theta_inc)+2*pi/period_x');	
var_Tm1.set('k_t_y_m1', '-sqrt(k0^2-k_t_x_m1^2)');	
var_Tm1.set('theta_t_m1', 'asin(sin(theta_inc)+lam0/period_x)');
var_Tm1.set('A_m1', 'intop2(acpr.p_t*exp(i*k_t_x_m1*x+i*k_t_y_m1*y))/period_x');	
var_Tm1.set('T_m1', '(abs(A_m1))^2*cos(theta_t_m1)/((abs(A_inc))^2*cos(theta_inc))');

var_Tm2 = comp1.variable.create('var_Tm2');
var_Tm2.label('Variables Tm2');
var_Tm2.set('k_t_x_m2', 'k0*sin(theta_inc)+4*pi/period_x');
var_Tm2.set('k_t_y_m2', '-sqrt(k0^2-k_t_x_m2^2)');
var_Tm2.set('theta_t_m2', 'asin(sin(theta_inc)+2*lam0/period_x)');
var_Tm2.set('A_m2', 'intop2(acpr.p_t*exp(i*k_t_x_m2*x+i*k_t_y_m2*y))/period_x');
var_Tm2.set('T_m2', '(abs(A_m2))^2*cos(theta_t_m2)/((abs(A_inc))^2*cos(theta_inc))');

var_R0 = comp1.variable.create('var_R0');
var_R0.label('Variables R0');
var_R0.set('k_r_x_0', 'k0*sin(theta_inc)');
var_R0.set('k_r_y_0', 'sqrt(k0^2-k_r_x_0^2)');
var_R0.set('theta_r_0', 'theta_inc');
var_R0.set('B_0', 'intop3(acpr.p_s*exp(i*k_r_x_0*x+i*k_r_y_0*y))/period_x');	
var_R0.set('R_0', '(abs(B_0))^2*cos(theta_r_0)/((abs(B_inc))^2*cos(theta_inc))');

var_R1 = comp1.variable.create('var_R1');
var_R1.label('Variables R1');
var_R1.set('k_r_x_1', 'k0*sin(theta_inc)-2*pi/period_x');
var_R1.set('k_r_y_1', 'sqrt(k0^2-k_r_x_1^2)');
var_R1.set('theta_r_1', 'asin(-sin(theta_inc)+lam0/period_x)');
var_R1.set('B_1', 'intop3(acpr.p_s*exp(i*k_r_x_1*x+i*k_r_y_1*y))/period_x');	
var_R1.set('R_1', '(abs(B_1))^2*cos(theta_r_1)/((abs(B_inc))^2*cos(theta_inc))');

var_R2 = comp1.variable.create('var_R2');
var_R2.label('Variables R2');
var_R2.set('k_r_x_2', 'k0*sin(theta_inc)-4*pi/period_x');
var_R2.set('k_r_y_2', 'sqrt(k0^2-k_r_x_2^2)');	
var_R2.set('theta_r_2', 'asin(-sin(theta_inc)+2*lam0/period_x)');
var_R2.set('B_2', 'intop3(acpr.p_s*exp(i*k_r_x_2*x+i*k_r_y_2*y))/period_x');
var_R2.set('R_2', '(abs(B_2))^2*cos(theta_r_2)/((abs(B_inc))^2*cos(theta_inc))');

var_Rm1 = comp1.variable.create('var_Rm1');
var_Rm1.label('Variables Rm1');
var_Rm1.set('k_r_x_m1', 'k0*sin(theta_inc)+2*pi/period_x');	
var_Rm1.set('k_r_y_m1', 'sqrt(k0^2-k_r_x_m1^2)');
var_Rm1.set('theta_r_m1', 'asin(-sin(theta_inc)-lam0/period_x)');
var_Rm1.set('B_m1', 'intop3(acpr.p_s*exp(i*k_r_x_m1*x+i*k_r_y_m1*y))/period_x');	
var_Rm1.set('R_m1', '(abs(B_m1))^2*cos(theta_r_m1)/((abs(B_inc))^2*cos(theta_inc))');

var_Rm2 = comp1.variable.create('var_Rm2');
var_Rm2.label('Variables Rm2');
var_Rm2.set('k_r_x_m2', 'k0*sin(theta_inc)+4*pi/period_x');
var_Rm2.set('k_r_y_m2', 'sqrt(k0^2-k_r_x_m2^2)');
var_Rm2.set('theta_r_m2', 'asin(-sin(theta_inc)-2*lam0/period_x)');
var_Rm2.set('B_m2', 'intop3(acpr.p_s*exp(i*k_r_x_m2*x+i*k_r_y_m2*y))/period_x');
var_Rm2.set('R_m2', '(abs(B_m2))^2*cos(theta_r_m2)/((abs(B_inc))^2*cos(theta_inc))');

var_tot = comp1.variable.create('var_tot');
var_tot.label('Variables Total');
var_tot.set('R_tot', 'R_0+R_1+R_2+R_m1+R_m2');
var_tot.set('T_tot', 'T_0+T_1+T_2+T_m1+T_m2');
var_tot.set('scat_tot', 'T_0+T_1+T_2+T_m1+T_m2+R_0+R_1+R_2+R_m1+R_m2');

% model.component('comp1').cpl('intop2').set('axisym', true);
% model.component('comp1').cpl('intop2').selection.geom('geom0', 1);
% model.component('comp1').cpl('intop2').selection.set([8]);

% % % 
% % % 
% mphsave(model,sprintf('%s/%i_model_ololo', dir_path.geom, sys_info.calc_tag),'component','on')
% 
% 1