function [optimumCandidate,optimumObjectiveValue,optimizationData] = optimization_differential_evolution(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,varargin)    
%OPTIMIZATION_DIFFERENTIAL_EVOLUTION technique stochastic optimization to
%   solve continuous and complex optimization problems. It uses mutation,
%   recombination and selection to optimize a given objective function.
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_DIFFERENTIAL_EVOLUTION(OBJECTIVEFUNCTION,INITIALDESIGN,
%	DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,CONSTRAINTFUNCTION) minimizes 
%	OBJECTIVEFUNCTION starting on INITIALDESIGN, on the design space defined by 
%	DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, subject to the constraint 
%	functions CONSTRAINTFUNCTION being less than or equal to zero, returning the
%	optimal design OPTIMUMCANDIDATE. 
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_DIFFERENTIAL_EVOLUTION(...NAME,VALUE,...) 
%   allows for the specification of additional options. These include:
%       - MaxIterations : maximum possible iterations
%       - PopulationSize : size of the Population
%       - MutationFactor : controls the scale of the difference between 
%       candidate solutions (individuals) when generating new solutions during 
%       the mutation step(0<F<1)
%       - CrossOverConstant : determines the probability of recombination between 
%       a donor vector and the current population vector to form a trial vector
%       (0<CR<1)
%       - 'SamplingFunction' : function handle on how to perform the sampling.
%		Default: @sampling_latin_hypercube.
%		- 'SamplingOptions' : options to be used with the sampling function.
%		Default is empty.
%       - 'ConvergenceCriterionFunction' : function handle for convergence
%       criteria
%       - 'ConvergenceCriterionOptions' : options to be used with
%       convergence function. Default is empty.
%		
%	[DESIGNOPTIMAL,OPTIMUMOBJECTIVEVALUE] = OPTIMIZATION_DIFFERENTIAL_EVOLUTION(...) 
%   also returns the value of the objective function OPTIMUMOBJECTIVEVALUE for 
%   the optimal design.
%
%	%[OPTIMUMCANDIDATE,OPTIMUMOBJECTIVEVALUE,OPTIMIZATIONDATA] =
%   OPTIMIZATION_DIFFERENTIAL_EVOLUTION(...) returns a structure with the 
%   processed parameters used during training. This can be later used for 
%   reproducibility, to check for any issues in the input, and/or for plotting 
%   the performance of the algorithm.  In particular:
%		- 'MutationRandomDonorIndex' : random individual index
%       - 'MutationInitialDesign' : initial point of Mutation for each design sample point
%       - 'MutationDirection' : direction of Mutation for each design sample point
%       - 'MutationDonorVector' : mutant vector for each design sample point
%       - 'RecombinationFactor' : random number generation for each design sample point
%       - 'RecombinationUseDonorVector' : donor vector indices for producing trial vector
%       - 'RecombinationTrialVector' : trial vector after recombination
%       - 'EvaluationTrialVectorObjectiveValue' : objective value of trial vector
%       - 'PopulationDesignNew' : updated population after selection
%       - 'PopulationObjectiveValueNew' : objective value of updated population
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
%		- 'MaxIterations' : integer
%		- 'PopulationSize' : integer
%		- 'MutationFactor' : double
%		- 'SamplingFunction' : function_handle
%		- 'SamplingOptions' : (1,nOptions) cell
%       - 'ConvergenceCriterionFunction' : function_handle
%       - 'ConvergenceCriterionOptions' : (1,nOptions) cell
%
%   Output:
%		- DESIGNOPTIMAL : (1,nDesignVariable) double
%		- OPTIMUMOBJECTIVEVALUE : double
%		- OPTIMIZATIONDATA : structure
%			-- MutationRandomDonorIndex : (PopulationSize,3) double
%			-- MutationInitialDesign : (PopulationSize,nDesignVariable) double
%			-- MutationDirection : (PopulationSize,nDesignVariable) double
%           -- MutationDonorVector : (PopulationSize,nDesignVariable) double
%           -- RecombinationFactor : (PopulationSize,1) double
%           -- RecombinationUseDonorVector : (PopulationSize,1) integer
%           -- RecombinationTrialVector : (PopulationSize,nDesignVariable) double
%           -- EvaluationTrialVectorObjectiveValue : (PopulationSize,1) double
%           -- PopulationDesignNew : (PopulationSize,nDesignVariable) double
%           -- PopulationObjectiveValueNew : (PopulationSize,1) double
%           -- OptimumCandidate : (1,nDesignVariable) double
%           -- OptimumCandidateObjectiveValue : (1,1) double
%           -- OptimumCandidateEvaluationData : (1,1) struct
%   
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Gaurav Vaibhav (Main Author)
%   Copyright 2025 Ali Abbas Kapadia (Contributor)
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

    % parsing inputs
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'PopulationSize', 40, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'MutationFactor', 0.8, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'CrossOverConstant', 0.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'PenaltyCoefficient', 1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p,'SamplingFunction',@sampling_random);
    addParameter(p,'SamplingOptions',{});
    addParameter(p,'ConvergenceCriterionFunction',@convergence_criterion_max_population_variance);
    addParameter(p,'ConvergenceCriterionOptions',{});
    parse(p, varargin{:});
    options = p.Results;

    nDimension = size(designSpaceLowerBound,2);
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];

    if isempty(@constraintFunction)
        constraintFunction = @(x) deal(zeros(size(x, 1), 1), zeros(size(x, 1), 1));
    end

    initialRng = rng;
    nSamplingGenerate = options.PopulationSize - size(initialDesign,1);

    population = [...
        initialDesign;...
        options.SamplingFunction(...
            designSpace,...
            nSamplingGenerate,...
            options.SamplingOptions{:});
        ];
     
    [objectiveValuePopulation,populationEvaluationData] = evaluate_optimization_objective(objectiveFunction,population,options.EvaluationObjectiveOptions{:});

    [optimumObjectiveValue ,iSelected] = min(objectiveValuePopulation);
    optimumCandidate = population(iSelected,:);
    optimumCandidateEvaluationData = populationEvaluationData(iSelected);

    % only log if necessary
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
            'PopulationEvaluationDataInitial',populationEvaluationData,...
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

        % Step 1: Mutation
        %DEL: you need to ensure that none of the entries are the same
        randomIndividualIndex = nan(options.PopulationSize,3);
        for i=1:options.PopulationSize
            randomIndividualIndex(i,:) = randperm(options.PopulationSize,3);
        end
        initialPoint = population(randomIndividualIndex(:,1),:);
        directionMutation = population(randomIndividualIndex(:,2),:) - population(randomIndividualIndex(:,3),:);
        maximumMutationFactor = region_limit_line_search([],initialPoint,directionMutation,designSpace);
        mutationFactor = min(maximumMutationFactor,options.MutationFactor);
        donorVector = initialPoint + mutationFactor.*directionMutation;

        % Step 2: Recombination
        trialRandomFactor = rand(options.PopulationSize,1);
        useDonorVector = (trialRandomFactor<=options.CrossOverConstant);

        trialVector = nan(options.PopulationSize,nDimension);
        trialVector(useDonorVector,:) = donorVector(useDonorVector,:);
        trialVector(~useDonorVector,:) = population(~useDonorVector,:);

        % Step 3: Evaluation
        % only evaluate designs which haven't been evaluated before
        trialVectorPointEvaluate = trialVector(useDonorVector,:);
        [trialVectorNewObjective,newEvaluationData] = evaluate_optimization_objective(objectiveFunction,trialVectorPointEvaluate,options.EvaluationObjectiveOptions{:});
        
        objectiveValueTrialVector = nan(options.PopulationSize,1);
        objectiveValueTrialVector(useDonorVector) = trialVectorNewObjective;
        objectiveValueTrialVector(~useDonorVector) = objectiveValuePopulation(~useDonorVector);

        % Step 4: Selection
        updatePopulation = (objectiveValueTrialVector<objectiveValuePopulation);     
        population(updatePopulation,:) = trialVector(updatePopulation,:);
        objectiveValuePopulation(updatePopulation,:) = objectiveValueTrialVector(updatePopulation,:);

        % check for the current optimum candidate
        [optimumObjectiveValue,iSelected] = min(objectiveValuePopulation);
        optimumCandidate = population(iSelected,:);
        if(useDonorVector(iSelected))
            iSelectEvaluation = convert_index_base(useDonorVector,iSelected,'forward');
            optimumCandidateEvaluationData = newEvaluationData(iSelectEvaluation);
        end
        
        % check for convergence based on information
        hasConverged = options.ConvergenceCriterionFunction(...
            optimumCandidate,optimumCandidatePrevious,...
            optimumObjectiveValue,optimumObjectiveValuePrevious,...
            population,populationPrevious,...
            objectiveValuePopulation,objectiveValuePopulationPrevious,...
            options.ConvergenceCriterionOptions{:});

        % log information if required
        if(isOutputOptimizationData)
            optimizationData.IterationData(iteration) = struct(...
                'MutationRandomDonorIndex',randomIndividualIndex,...
                'MutationInitialDesign',initialPoint,...
                'MutationDirection',directionMutation,...
                'MutationDonorVector',donorVector,...
                'RecombinationFactor',trialRandomFactor,...
                'RecombinationUseDonorVector',useDonorVector,...
                'RecombinationTrialVector',trialVector,...
                'EvaluationTrialVectorObjectiveValue',objectiveValueTrialVector,...
                'EvaluationNewEvaluationData',newEvaluationData,...
                'PopulationDesignNew', population,...
                'PopulationObjectiveValueNew',objectiveValuePopulation,...
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

