function params = LoadParameters(params)

% Incident wave
params.wave.p0 = 1;
params.wave.theta_inc = 0;
params.wave.freq0 = 5 * 1e6;
params.wave.c_water = 1500;
params.wave.lam0 = params.wave.c_water/params.wave.freq0;
params.wave.ff_zone = 6*params.wave.lam0;


% Metalens
params.metalens.focus = 5*1e-3;
params.metalens.x0 = 0.3*1e-3;
params.metalens.N = 1;

% Unit cell size
[p, ~] = GratingPeriod(params);
params.unit_cell.period_x = p(params.metalens.N);
params.unit_cell.min_element = params.unit_cell.period_x/10;
params.unit_cell.ax = params.unit_cell.period_x - 2*params.unit_cell.min_element;
% params.unit_cell.ay = 1.5*params.unit_cell.period_x;

% Host 
params.host.size_x = params.unit_cell.period_x;

% Elements
params.elements.n_vertices = 4;
params.elements.vx = {};
params.elements.vy = {};


% PML
params.pml.width = params.wave.ff_zone/2;
params.pml.scaling = 1;
params.pml.curvature = 1;



% Mesh
params.mesh.max_element = params.wave.lam0/6;
params.mesh.min_element = params.wave.lam0/100;
params.mesh.narrow_resolution = 5;
params.mesh.curvature = 0.3;
params.mesh.growth_rate = 1.5;
params.mesh.pml_distribution = 10;
params.mesh.save_mesh = false;


%% Algorithm parameters
params.algorithm.n_population = 6;
params.algorithm.max_iter = 30;
params.algorithm.n_offsprings = 1;
params.algorithm.mutation_rate = 0.25;
params.algotithm.beta = 0.1;
params.algorithm.sigma0 = 0.1;
params.algorithm.gamma0 = 0.1;

