function output = RealGeneticAlgorithm(problem, params, sys_info)

CostFun = problem.CostFun;
n_var = problem.n_var;
VarMin = problem.VarMin;
VarMax = problem.VarMax;
n_population = params.algorithm.n_population;
max_iter = params.algorithm.max_iter;
n_offsprings = params.algorithm.n_offsprings;
mutation_rate = params.algorithm.mutation_rate;
sigma0 = params.algorithm.sigma0;
gamma0 = params.algorithm.gamma0;
n_children = round(n_population*n_offsprings/2) * 2;

% Initialization
individ.position = [];
individ.cost = [];

best_sol.cost = inf;

sys_info.calc_tag = 0;
population = repmat(individ, n_population, 1);
for idx_pop = 1:n_population
    sys_info.calc_tag = sys_info.calc_tag + 1;
    population(idx_pop).position = unifrnd(VarMin, VarMax, [1, n_var]);
%     population(idx_pop).cost = CostFun(population(idx_pop).position);
    population(idx_pop).cost = ...
        metagrating(population(idx_pop).position, params, sys_info);
    if population(idx_pop).cost < best_sol.cost
        best_sol = population(idx_pop);
    end
end

% Store the best cost
best_cost = nan(max_iter, 1);
best_tag = 0;

figure();

% Evolution loop
for idx_iter = 1:max_iter

    population_offspring = repmat(individ, n_children/2, 2);
    
    % Selection probabilities
    c = [population.cost];
    avgc = mean(c);
    if avgc ~= 0
        c = c/avgc;
    end
    probs = exp(-params.algotithm.beta * c);

    % Crossover
    for idx_k = 1:n_children/2

%         % Random selection of parents
%         q = randperm(n_population);
%         parent_1 = population(q(1));
%         parent_2 = population(q(2));
        parent_1 = population(RouletteWheelSelection(probs));
        parent_2 = population(RouletteWheelSelection(probs));
        
        % Crossover
        [population_offspring(idx_k,1).position, ...
            population_offspring(idx_k,2).position] = ...
            UniformCrossover(parent_1.position, parent_2.position, gamma0);
    end

    population_offspring = population_offspring(:);

    % Mutation
    for l = 1:n_children

        sys_info.calc_tag = sys_info.calc_tag + 1;

        population_offspring(l).position = ...
            Mutate(population_offspring(l).position, mutation_rate, sigma0);
        vx_array = population_offspring(l).position(1:2:end);
        vy_array = population_offspring(l).position(2:2:end);
        idx_error = 1;
        dx = cell(1,length(params.elements.n_vertices));
        while idx_error ~= 0
            idx_error = 0;
            for idx_v = 1:length(params.elements.n_vertices)
                dx{idx_v} = sqrt((vx_array - vx_array(idx_v)).^2 + (vy_array - vy_array(idx_v)).^2);

                dx_x{idx_v} = abs(vx_array - vx_array(idx_v));
                dx_y{idx_v} = abs(vy_array - vy_array(idx_v));
                
                if dx{idx_v} < params.unit_cell.min_element
                    idx_error = idx_error + 1;
                    vx_array(idx_v) = vx_array(idx_v) + params.unit_cell.min_element;
                    vy_array(idx_v) = vy_array(idx_v) + params.unit_cell.min_element;
                elseif dx_x{idx_v} < params.unit_cell.min_element
                    idx_error = idx_error + 1;
                    vx_array(idx_v) = vx_array(idx_v) + params.unit_cell.min_element;
                elseif dx_y{idx_v} < params.unit_cell.min_element
                    idx_error = idx_error + 1;
                    vy_array(idx_v) = vy_array(idx_v) + params.unit_cell.min_element;
                end
            end
        end
        population_offspring(l).position(1:2:end) = vx_array;
        population_offspring(l).position(2:2:end) = vy_array;
        population_offspring(l).position = ...
            max(population_offspring(l).position, VarMin);
        population_offspring(l).position = ...
            min(population_offspring(l).position, VarMax);
%         population_offspring(l).cost = ...
%             CostFun(population_offspring(l).position);
        population_offspring(l).cost = ...
            metagrating(population_offspring(l).position, params, sys_info);
        if population_offspring(l).cost < best_sol.cost
            best_sol = population_offspring(l);
            best_tag = sys_info.calc_tag;
        end
    end
    best_cost(idx_iter) = best_sol.cost;

    semilogy(1:max_iter, best_cost, 'k', 'linewidth', 2);
    title_str = ...
        sprintf('Best cost = %f, design %i', best_cost(idx_iter), best_tag);
    title(title_str);
    xlabel('Iteration')
    ylabel('Cost')
    drawnow;

    % Merging
    population = [population; population_offspring];

    % Sorting
    [~, sort_order] = sort([population.cost]);
    population = population(sort_order);


    % Selection
    population = population(1:n_population);
    
    fclose('all');
end

output.population = population;
output.best_sol = best_sol;
output.best_cost = best_cost;
