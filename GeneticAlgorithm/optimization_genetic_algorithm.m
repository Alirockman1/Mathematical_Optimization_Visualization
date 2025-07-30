function [optimumCandidate,optimumObjectiveValue,optimizationData] = optimization_genetic_algorithm(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,varargin)
%OPTIMIZATION_GENETIC_ALGORITHM technique stochastic optimization to
%   solve continuous and complex optimization problems. It uses mutation,
%   recombination and selection to optimize a given objective function.
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_GENETIC_ALGORITHM(OBJECTIVEFUNCTION,INITIALDESIGN,
%	DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,CONSTRAINTFUNCTION) minimizes 
%	OBJECTIVEFUNCTION starting on INITIALDESIGN, on the design space defined by 
%	DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, subject to the constraint 
%	functions CONSTRAINTFUNCTION being less than or equal to zero, returning the
%	optimal design OPTIMUMCANDIDATE. 
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_GENETIC_ALGORITHM(...NAME,VALUE,...) 
%   allows for the specification of additional options. These include:
%       - MaxIterations : maximum possible iterations
%       - PopulationSize : size of the Population
%       - FitnessFunction : function handle of evaluation
%       - FitnessOptions : options to be used with the fitness function
%       - SelectionFunction : function handle of selection
%       - SelectionOptions : options to be used with the selection function
%       - CrossoverRate : probability of cross over
%       - CrossoverFunction : function handle of cross over
%       - CrossoverOptions : options to be used with the crossover function
%       - MutationRate : probability of mutation
%       - MutationFunction : function handle of mutation
%       - MutationOptions : options to be used with the mutation function
%       - 'SamplingFunction' : function handle on how to perform the sampling.
%		Default: @sampling_latin_hypercube.
%		- 'SamplingOptions' : options to be used with the sampling function.
%		Default is empty.
%       - 'ConvergenceCriterionFunction' : function handle for convergence
%       criteria
%       - 'ConvergenceCriterionOptions' : options to be used with
%       convergence function. Default is empty.
%		
%	[DESIGNOPTIMAL,OPTIMUMOBJECTIVEVALUE] = OPTIMIZATION_GENETIC_ALGORITHM(...) 
%   also returns the value of the objective function OPTIMUMOBJECTIVEVALUE for 
%   the optimal design.
%
%	%[OPTIMUMCANDIDATE,OPTIMUMOBJECTIVEVALUE,OPTIMIZATIONDATA] =
%   OPTIMIZATION_GENETIC_ALGORITHM(...) returns a structure with the 
%   processed parameters used during training. This can be later used for 
%   reproducibility, to check for any issues in the input, and/or for plotting 
%   the performance of the algorithm.  In particular:
%		- 'FitnessValue' : objective function value
%       - 'SelectionIndex' : indices of population for selection according to selection function
%       - 'SelectionPopulation' : updated population after selection
%       - 'CrossoverPopulation' : population formed after cross over and random number comparison
%       - 'MutationFactor' : scaling factor that determines how much change is applied to the mutated values
%       - 'MutationIsMutated' : population indices giving mutation occurence
%       - 'MutationPopulation' : population formed after mutation
%       - 'EvaluationTrialVectorObjectiveValue' : objective value of trial vector
%       - 'EvaluationObjectiveValue' : objective value of updated population
%       - 'OptimumCandidate' : optimal design variable for current iteration 
%       - 'OptimumCandidateObjectiveValue' : optimal objective value for current iteration
%       - 'OptimumCandidateEvaluationData' : optimal values for design, point, base and 
%       unconstrained objective values, equality and inequality constraint penalty and barrier values
%
%   Input:
%		- OBJECTIVEFUNCTION : function_handle
%		- INITIALDESIGN : (1,nDesignVariable) double
%		- DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%		- DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%		- CONSTRAINTFUNCTION : function_handle
%		- 'MaxIterations' : integer
%		- 'PopulationSize' : integer
%		- 'FitnessFunction' : function_handle
%       - 'FitnessOptions' : (1,nOptions) cell
%		- 'SelectionFunction' : function_handle
%       - 'SelectionOptions' : (1,nOptions) cell
%       - 'CrossoverRate' : double
%		- 'CrossoverFunction' : function_handle
%       - 'CrossoverOptions' : (1,nOptions) cell
%		- 'MutationRate' : double
%		- 'MutationFunction' : function_handle
%       - 'MutationOptions' : (1,nOptions) cell
%		- 'SamplingFunction' : function_handle
%		- 'SamplingOptions' : (1,nOptions) cell
%       - 'ConvergenceCriterionFunction' : function_handle
%       - 'ConvergenceCriterionOptions' : (1,nOptions) cell
%
%   Output:
%		- DESIGNOPTIMAL : (1,nDesignVariable) double
%		- OPTIMUMOBJECTIVEVALUE : double
%		- OPTIMIZATIONDATA : structure
%			-- FitnessValue : (PopulationSize,1) double
%			-- SelectionIndex : (PopulationSize,1) double
%			-- SelectionPopulation : (PopulationSize,nDesignVariable) double
%           -- CrossoverPopulation : (PopulationSize,nDesignVariable) double
%           -- MutationFactor : (PopulationSize,1) double
%           -- MutationIsMutated : (PopulationSize,1) boolean
%           -- MutationPopulation : (PopulationSize,nDesignVariable) double
%           -- EvaluationObjectiveValue : (PopulationSize,1) double
%           -- OptimumCandidate : (1,nDesignVariable) double
%           -- OptimumCandidateObjectiveValue : (1,1) double
%           -- OptimumCandidateEvaluationData : (1,1) struct
%   
%   Copyright 2024 Gaurav Vaibhav, Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%   
%       http://www.apache.org/licenses/LICENSE-2.0
%   
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.    

    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'PopulationSize', 40, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'FitnessFunction',@genetic_algorithm_fitness_inverse);
    addParameter(p, 'FitnessOptions', {});
    addParameter(p, 'SelectionFunction',@genetic_algorithm_selection_roulette_wheel);
    addParameter(p, 'SelectionOptions', {});
    addParameter(p, 'CrossoverRate', 0.92);
    addParameter(p, 'CrossoverFunction', @genetic_algorithm_crossover_variable_exchange);
    addParameter(p, 'CrossoverOptions', {});
    addParameter(p, 'MutationRate', 0.2);
    addParameter(p, 'MutationFunction', @genetic_algorithm_mutation_uniform);
    addParameter(p, 'MutationOptions', {});
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'SamplingFunction',@sampling_random);
    addParameter(p, 'SamplingOptions',{});
    addParameter(p, 'ConvergenceCriterionFunction',@convergence_criterion_max_population_variance);
    addParameter(p, 'ConvergenceCriterionOptions',{});
    parse(p, varargin{:});
    options = p.Results;

    nDimension = size(designSpaceLowerBound,2);
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];

    nSamplingGenerate = options.PopulationSize - size(initialDesign,1);

    % generating initial population
    initialRng = rng;
    population = [...
        initialDesign;...
        options.SamplingFunction(...
            designSpace,...
            nSamplingGenerate,...
            options.SamplingOptions{:});
        ];
    
%     objectiveValuePopulation = objectiveFunction(population);
    [objectiveValuePopulation,populationEvaluationData] = evaluate_optimization_objective(objectiveFunction,population,options.EvaluationObjectiveOptions{:});

    [optimumObjectiveValue ,iSelected] = min(objectiveValuePopulation);
    optimumCandidate = population(iSelected,:);
    optimumCandidateEvaluationData = populationEvaluationData(iSelected);

    %DEL: only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
        optimizationData.ProblemData = struct(...
            'ObjectiveFunction',objectiveFunction,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options,...
            'InitialRngState',initialRng);

        optimizationData.InitialData = struct(...
            'PopulationDesignInitial', population,...
            'PopulationObjectiveValueInitial',objectiveValuePopulation,...
            'OptimumCandidate',optimumCandidate,...
            'OptimumCandidateObjectiveValue',optimumObjectiveValue,...
            'OptimumCandidateEvaluationData',optimumCandidateEvaluationData);
    end
    
    hasConverged = false;

    % adding for easy future reference
    populationPrevious = population;
    objectiveValuePopulationPrevious = objectiveValuePopulation;
    optimumCandidatePrevious = optimumCandidate;
    optimumObjectiveValuePrevious = optimumObjectiveValue;

    iteration = 1;
    while(~hasConverged && iteration <= options.MaxIterations) 
        % Step 1: Fitness
        fitnessValue = options.FitnessFunction(objectiveValuePopulation,options.FitnessOptions{:}); 

        % Step 2: Selection
        selectionIndex = options.SelectionFunction(fitnessValue,options.SelectionOptions{:});
        selectedPopulation = population(selectionIndex, :);

        % Step 2: Crossover
        offspringPopulation = selectedPopulation;
        for i = 1:2:(options.PopulationSize-1)
            if rand < options.CrossoverRate
                offspringPopulation([i;i+1], :) = options.CrossoverFunction(...
                    selectedPopulation([i;i+1],:),fitnessValue([i;i+1]),options.CrossoverOptions{:});
            end
        end

        % Step 3: Mutation
        mutationFactor = rand(options.PopulationSize,1);
        isMutate = (mutationFactor<options.MutationRate);
        mutatedPopulation = offspringPopulation;
        mutatedPopulation(isMutate,:) = options.MutationFunction(offspringPopulation(isMutate,:),designSpace,options.MutationOptions{:});

        % Step 4: Evaluation
        population = mutatedPopulation;
        [objectiveValuePopulation,newEvaluationData] = evaluate_optimization_objective(objectiveFunction,population,options.EvaluationObjectiveOptions{:});
        
        % check for the current optimum candidate
        [optimumObjectiveValue ,iSelected] = min(objectiveValuePopulation);
        optimumCandidate = population(iSelected,:);
        optimumCandidateEvaluationData = newEvaluationData(iSelected);

        hasConverged = options.ConvergenceCriterionFunction(...
            optimumCandidate,optimumCandidatePrevious,...
            optimumObjectiveValue,optimumObjectiveValuePrevious,...
            population,populationPrevious,...
            objectiveValuePopulation,objectiveValuePopulationPrevious,...
            options.ConvergenceCriterionOptions{:});

        % log information if required
        if(isOutputOptimizationData)
            optimizationData.IterationData(iteration) = struct(...
                'FitnessValue',fitnessValue,...
                'SelectionIndex',selectionIndex,...
                'SelectionPopulation',selectedPopulation,...
                'CrossoverPopulation',offspringPopulation,...
                'MutationFactor',mutationFactor,...
                'MutationIsMutated',isMutate,...
                'MutationPopulation',mutatedPopulation,...
                'PopulationDesignNew',population,...
                'EvaluationObjectiveValue',objectiveValuePopulation,...
                'OptimumCandidate',optimumCandidate,...
                'OptimumCandidateObjectiveValue',optimumObjectiveValue,...
                'OptimumCandidateEvaluationData',optimumCandidateEvaluationData);
        end

        % update state
        populationPrevious = population;
        objectiveValuePopulationPrevious = objectiveValuePopulation;
        optimumCandidatePrevious = optimumCandidate;
        optimumObjectiveValuePrevious = optimumObjectiveValue;
        iteration = iteration + 1;
    end
end