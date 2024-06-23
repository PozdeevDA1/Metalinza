function output = BinaryGeneticAlgorithm(problem, params, sys_info)

CostFun = problem.CostFun;
n_var = problem.n_var;
n_population = params.algorithm.n_population;
max_iter = params.algorithm.max_iter;
n_offsprings = params.algorithm.n_offsprings;
mutation_rate = params.algorithm.mutation_rate;


n_children = round(n_population*n_offsprings/2) * 2;

% Initialization
individ.position = [];
individ.cost = [];

best_sol.cost = inf;

sys_info.calc_tag = 0;
population = repmat(individ, n_population, 1);
% population(1).position = zeros(1, n_var);
% population(1).cost = ...
%         defects_wall_GA(population(1).position, params, sys_info);
for idx_pop = 1:n_population
    sys_info.calc_tag = sys_info.calc_tag + 1;
    population(idx_pop).position = randi([0, 1], 1, n_var);


%     population(idx_pop).cost = CostFun(population(idx_pop).position);
    population(idx_pop).cost = ...
        metagrating(population(idx_pop).position, params, sys_info);
    if population(idx_pop).cost < best_sol.cost
        best_sol = population(idx_pop);
    end
end

% Store the best cost
best_cost = nan(max_iter, 1);

% Evolution loop
for idx_iter = 1:max_iter
    disp(idx_iter)
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
            Crossover(parent_1.position, parent_2.position);
    end

    population_offspring = population_offspring(:);

    % Mutation
    for l = 1:n_children

        sys_info.calc_tag = sys_info.calc_tag + 1;

        population_offspring(l).position = ...
            Mutate(population_offspring(l).position, mutation_rate);
%         population_offspring(l).cost = ...
%             CostFun(population_offspring(l).position);
        population_offspring(l).cost = ...
            metagrating(population_offspring(l).position, params, sys_info);
        if population_offspring(l).cost < best_sol.cost
            best_sol = population_offspring(l);
        end
    end
    best_cost(idx_iter) = best_sol.cost;
    
    if idx_iter > 50
       if  best_cost(idx_iter) == best_cost(idx_iter-50)
           break;
       end
    end
    % Merging
    population = [population; population_offspring];

    % Sorting
    [~, sort_order] = sort([population.cost]);
    population = population(sort_order);


    % Selection
    population = population(1:n_population);
end

output.population = population;
output.best_sol = best_sol;
output.best_cost = best_cost;
