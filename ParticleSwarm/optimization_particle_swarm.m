function[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_particle_swarm(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,varargin)
%OPTIMIZATION_PARTICLE_SWARM Perform optimization using Particle Swarm Optimization (PSO).
%
%   General Description:
%   This function applies the Particle Swarm Optimization (PSO) algorithm to minimize a given
%   objective function over a specified design space. Inspired by the social behavior of birds,
%   PSO explores the design space by iteratively updating candidate solutions based on personal
%   and global best knowledge, with optional convergence checks and parameter customization.
%
%   Inputs:
%       - objectiveFunction : Function handle to the objective function.
%       - initialDesign : Initial population seed, (nInitial, nVariables) double.
%       - designSpaceLowerBound : Lower bounds of the design variables, (1, nVariables) double.
%       - designSpaceUpperBound : Upper bounds of the design variables, (1, nVariables) double.
%       - Name-Value Pairs:
%           - 'MaxIterations' : Maximum number of iterations (default: 100).
%           - 'PopulationSize' : Total number of particles (default: 10).
%           - 'InertiaCoefficient' : Momentum of particle movement (default: 0.5).
%           - 'CognitiveCoefficient' : Influence of personal best (default: 0.5).
%           - 'SocialCoefficient' : Influence of global best (default: 0.5).
%           - 'VelocityBounds' : Velocity limits, scalar or [2 x nVariables] (default: 0.3).
%           - 'UseUniformlyDistributedCoefficients' : Boolean for dynamic coefficient scaling (default: false).
%           - 'SamplingFunction' : Function handle for sampling (default: @sampling_random).
%           - 'SamplingOptions' : Cell array of arguments for sampling function (default: {}).
%           - 'EvaluationObjectiveOptions' : Options passed to the objective evaluation (default: {}).
%           - 'ConvergenceCriterionFunction' : Handle for convergence function (default: @convergence_criterion_max_population_variance).
%           - 'ConvergenceCriterionOptions' : Cell array of arguments for convergence function (default: {}).
%
%   Outputs:
%       - optimumCandidate : Best design found, (1, nVariables) double.
%       - optimumObjectiveValue : Objective value at the optimumCandidate.
%       - optimizationData : Struct containing full algorithm data:
%           - ProblemData : Initial setup and configuration.
%           - InitialData : Starting design, velocity, and initial evaluations.
%           - IterationData : Struct array per iteration with fields:
%               * Velocities
%               * PopulationDesignNew
%               * PopulationObjectiveValue
%               * PopulationWasImproved
%               * IndividualBestPosition
%               * IndividualBestObjectiveValue
%               * OptimumCandidate
%               * OptimumCandidateObjectiveValue
%               * OptimumCandidateEvaluationData
%
%   Functionality:
%       The algorithm initializes a swarm of particles and iteratively updates their
%       positions and velocities using inertia, cognitive, and social influences.
%       The best solution is tracked across all particles and iterations.
%       The function supports early stopping based on a convergence criterion.
%
%   Dependencies:
%       - evaluate_optimization_objective
%       - region_limit_line_search
%       - sampling_random or custom sampling function
%       - convergence_criterion_max_population_variance or custom convergence function
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
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

    %Parsing inputs
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'PopulationSize', 10, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'InertiaCoefficient', 0.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'CognitiveCoefficient', 0.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'SocialCoefficient', 0.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'VelocityBounds', 0.3, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'UseUniformlyDistributedCoefficients',false);
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'SamplingFunction',@sampling_random);
    addParameter(p, 'SamplingOptions',{});
    addParameter(p, 'ConvergenceCriterionFunction',@convergence_criterion_max_population_variance);
    addParameter(p, 'ConvergenceCriterionOptions',{});

    parse(p, varargin{:});
    options = p.Results;
    
    nDimension = size(designSpaceLowerBound,2);
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];

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
    individualBestPosition = population;
    individualBestObjectiveValue = objectiveValuePopulation;

    [optimumObjectiveValue ,iSelected] = min(individualBestObjectiveValue);
    optimumCandidate = individualBestPosition(iSelected,:);
    optimumCandidateEvaluationData = populationEvaluationData(iSelected);

    if(isscalar(options.VelocityBounds))
        options.VelocityBounds = abs(options.VelocityBounds).*[-ones(1,nDimension);ones(1,nDimension)];
    elseif(size(options.VelocityBounds,1)==1)
        options.VelocityBounds = [-abs(options.VelocityBounds);abs(options.VelocityBounds)];
    end

    velocities = options.SamplingFunction(...
        options.VelocityBounds,...
        options.PopulationSize,...
        options.SamplingOptions{:});

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
            'OptimumCandidate',optimumCandidate,...
            'OptimumCandidateObjectiveValue',optimumObjectiveValue,...
            'OptimumCandidateEvaluationData',optimumCandidateEvaluationData,...
            'Velocities',velocities);
    end
   
    hasConverged = false;

    % adding for easy future reference
    populationPrevious = population;
    objectiveValuePopulationPrevious = objectiveValuePopulation;
    optimumCandidatePrevious = optimumCandidate;
    optimumObjectiveValuePrevious = optimumObjectiveValue;

    iteration = 1;
    while(~hasConverged && iteration <= options.MaxIterations)
        if(options.UseUniformlyDistributedCoefficients)
            bilateralUniformDistribution = @(x) -x + rand(options.PopulationSize,1).*(2*x);
            cognitiveFactor = bilateralUniformDistribution(options.CognitiveCoefficient);
            socialFactor = bilateralUniformDistribution(options.SocialCoefficient);
        else
            cognitiveFactor = options.CognitiveCoefficient;
            socialFactor = options.SocialCoefficient;
        end

        % Update velocities
        velocities = ...
            options.InertiaCoefficient.*velocities + ...
            cognitiveFactor.*(individualBestPosition-population) + ...
            socialFactor.*(optimumCandidate-population);
        velocities = max(min(velocities,options.VelocityBounds(2,:)), options.VelocityBounds(1,:));

        % velocity limitation - don't go outside design space
        velocitiesFactorMax = region_limit_line_search([],population,velocities,designSpace);
        velocities = velocities.*min(velocitiesFactorMax,1);

        % update population
        population = population + velocities;

        % evaluate
        [objectiveValuePopulation,newEvaluationData] = evaluate_optimization_objective(objectiveFunction,population,options.EvaluationObjectiveOptions{:});

        % find where there was improvement
        isImproved = (objectiveValuePopulation<individualBestObjectiveValue);
        individualBestObjectiveValue(isImproved,:) = objectiveValuePopulation(isImproved,:);
        individualBestPosition(isImproved,:) = population(isImproved,:);

        % update overall information
        [optimumObjectiveValue ,iSelected] = min(individualBestObjectiveValue);
        optimumCandidate = individualBestPosition(iSelected,:);
        optimumCandidateEvaluationData = newEvaluationData(iSelected);

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
                'Velocities',velocities,...
                'PopulationDesignNew',population,...
                'PopulationObjectiveValue',objectiveValuePopulation,...
                'PopulationWasImproved',isImproved,...
                'IndividualBestPosition',individualBestPosition,...
                'IndividualBestObjectiveValue',individualBestObjectiveValue,...
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
        