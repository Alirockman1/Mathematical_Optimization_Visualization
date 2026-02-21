function [optimumCandidate,optimumObjectiveValue,optimizationData] = optimization_barrier_penalty(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,constraints,varargin)
%OPTIMIZATION_BARRIER_PENALTY technique to solve complex constrained optimization
%   problems with stochastic optimization methods by introducing penalty. 
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_BARRIER_PENALTY(OBJECTIVEFUNCTION,INITIALDESIGN,
%	DESIGNSPACELOWERBOUND,DESIGNSPACEUPPERBOUND,CONSTRAINTS) minimizes 
%	OBJECTIVEFUNCTION starting on INITIALDESIGN, on the design space defined by 
%	DESIGNSPACELOWERBOUND and DESIGNSPACEUPPERBOUND, subject to the constraint 
%	functions CONSTRAINTS being less than or equal to zero, returning the
%	optimal design OPTIMUMCANDIDATE. 
%
%	OPTIMUMCANDIDATE = OPTIMIZATION_BARRIER_PENALTY(...NAME,VALUE,...) 
%   allows for the specification of additional options. These include:
%       - InitialPenaltyFactor : factor for initial penalty
%       - PenaltyFactorUpdateFunction : function for updating penalty factor
%       - PenaltyFactorUpdateOptions : function for updating penalty factor
%       - MaxIterations : maximum possible iterations
%       - BaseOptimizationFunction : stochastic optimization technique
%       - BaseOptimizationOptions : options associated with above optimization 
%       technique
%       - 'ConvergenceCriterionFunction' : function handle for convergence criteria
%       - 'ConvergenceCriterionOptions' : options to be used with convergence 
%       function. Default is empty.
%		
%	[OPTIMUMCANDIDATE,OPTIMUMOBJECTIVEVALUE] = OPTIMIZATION_DIFFERENTIAL_EVOLUTION(...) 
%   also returns the value of the objective function OPTIMUMOBJECTIVEVALUE for 
%   the optimal design.
%
%	%[OPTIMUMCANDIDATE,OPTIMUMOBJECTIVEVALUE,OPTIMIZATIONDATA] =
%   OPTIMIZATION_DIFFERENTIAL_EVOLUTION(...) returns a structure OPTIMIZATIONDATA 
%   with the processed parameters used during training. This can be later used for 
%   reproducibility, to check for any issues in the input, and/or for plotting 
%   the performance of the algorithm.  In particular:
%       - 'OptimumCandidate' : optimal design variable for current iteration 
%       - 'OptimumCandidateObjectiveValue' : optimal objective value for current iteration
%       - 'BaseOptimizationData' : output optimization data from optimization function 
%       containing progress of variables in each iteration 
%       - 'PenaltyFactor' : penalty factor for current iteration 
% 
%   Input:
%		- OBJECTIVEFUNCTION : function_handle
%		- INITIALDESIGN : (1,nDesignVariable) double
%		- DESIGNSPACELOWERBOUND : (1,nDesignVariable) double
%		- DESIGNSPACEUPPERBOUND : (1,nDesignVariable) double
%		- CONSTRAINTS : function_handle
%		- InitialPenaltyFactor' : double
%		- 'PenaltyFactorUpdateFunction' : function_handle
%		- 'PenaltyFactorUpdateOptions' : (1,nOptions) cell
%		- 'MaxIterations' : integer
%		- 'BaseOptimizationFunction' : function_handle
%		- 'BaseOptimizationOptions' : (1,nOptions) cell
%       - 'ConvergenceCriterionFunction' : function_handle
%       - 'ConvergenceCriterionOptions' : (1,nOptions) cell
%
%   Output:
%		- OPTIMUMCANDIDATE : (1,nDesignVariable) double
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
%       - PENALTYFACTOR : (iteration,1) double
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
    addParameter(p,'InitialPenaltyFactor', 1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'PenaltyFactorUpdateFunction', @penalty_factor_update);
    addParameter(p,'PenaltyFactorUpdateOptions', {});
    addParameter(p,'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'BaseOptimizationFunction',@optimization_differential_evolution);
    addParameter(p,'BaseOptimizationOptions',{});
    addParameter(p,'ConvergenceCriterionFunction',@convergence_criterion_optimum_candidate_variance);
    addParameter(p,'ConvergenceCriterionOptions',{});
    parse(p, varargin{:});
    options = p.Results;
    
    initialPenalty = options.InitialPenaltyFactor;
    penaltyFactor = initialPenalty;

    optimumCandidatePrevious = initialDesign;
    optimumObjectiveValuePrevious = +inf;

    % only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
        optimizationData.ProblemData = struct(...
            'ObjectiveFunction',objectiveFunction,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options,...
            'InitialRngState',rng);

        optimizationData.InitialData = struct(...
            'DesignInitial', initialDesign,...
            'ObjectiveValueInitial',+inf,...
            'OptimumCandidate',initialDesign,...
            'OptimumCandidateObjectiveValue',+inf);
    end

    iteration = 1;
    hasConverged = false;
	while(~hasConverged && iteration<options.MaxIterations)
        evaluationObjectiveOptions = {'ConstraintFunctions',constraints,'PenaltyConstant',penaltyFactor};
        [~,baseOptimizationOptions] = merge_name_value_pair_argument(options.BaseOptimizationOptions,{'EvaluationObjectiveOptions',evaluationObjectiveOptions});

		[optimumCandidate,optimumObjectiveValue,baseOptimizationData] = options.BaseOptimizationFunction(objectiveFunction,initialDesign,designSpaceLowerBound,designSpaceUpperBound,...
            baseOptimizationOptions{:});

   		% Updating Penalty Factor
        penaltyFactor = options.PenaltyFactorUpdateFunction(initialPenalty,penaltyFactor,iteration,options.PenaltyFactorUpdateOptions{:});

		% logging the results
        optimizationData.IterationData(iteration) = struct(...
			'OptimumCandidate',optimumCandidate,...
			'OptimumCandidateObjectiveValue',optimumObjectiveValue,...
			'BaseOptimizationData',baseOptimizationData,...
			'PenaltyFactor',penaltyFactor);
        iteration = iteration + 1;

         % check for convergence based on information
        hasConverged = options.ConvergenceCriterionFunction(...
            optimumCandidate,optimumCandidatePrevious,...
            optimumObjectiveValue,optimumObjectiveValuePrevious,...
            optimumCandidate,optimumCandidatePrevious,...
            optimumObjectiveValue,optimumObjectiveValuePrevious,...
            options.ConvergenceCriterionOptions{:});
      
        % update state
        optimumObjectiveValuePrevious = optimumObjectiveValue;
        optimumCandidatePrevious = optimumCandidate;
        initialDesign = optimumCandidate;
	end
end