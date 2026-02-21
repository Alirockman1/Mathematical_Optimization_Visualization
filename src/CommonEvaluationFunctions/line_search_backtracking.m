function [stepSize, stepSizeData, lineSearchData] = line_search_backtracking(objectiveFunction, design, searchDirection, varargin)
%OPTIMIZATION_STEEPEST_DESCENT Performs gradient-based optimization using steepest descent with line search.
%
%   Syntax:
%       optimumCandidate = optimization_steepest_descent(objectiveFunction, initialDesign, lowerBound, upperBound)
%       [optimumCandidate, optimumObjectiveValue] = ...
%       [optimumCandidate, optimumObjectiveValue, optimizationData] = ...
%
%   Description:
%       This function implements a constrained steepest descent optimization method.
%       It iteratively updates a design vector using the negative gradient direction and
%       a line search step size. Gradient evaluation and objective evaluations are modular
%       and support additional runtime configuration.
%
%   Inputs:
%       objectiveFunction           : Function handle to the objective function.
%       initialDesign               : Initial design vector (1 x n or n x 1).
%       designSpaceLowerBound       : Lower bounds for each design variable (1 x n).
%       designSpaceUpperBound       : Upper bounds for each design variable (1 x n).
%
%   Name-Value Pair Arguments (optional):
%       'MaxIterations'             : Maximum number of iterations (default: 100).
%       'Tolerance'                 : Convergence threshold (default: 1e-6).
%       'UseNormalizedGradientLineSearch' : Normalize search direction (default: false).
%       'LineSearchFunction'        : Handle to line search function (default: @line_search_backtracking).
%       'EvaluationObjectiveOptions': Cell array of options for objective evaluation.
%       'GradientOptions'           : Cell array of options for gradient computation.
%       'LineSearchOptions'         : Cell array of options for line search.
%       'ConvergenceCriterionFunction' : Custom convergence function (default: @convergence_criterion_optimum_candidate_variance).
%       'ConvergenceCriterionOptions'  : Cell array of options for convergence criterion.
%
%   Outputs:
%       optimumCandidate            : Final design vector that minimizes the objective.
%       optimumObjectiveValue       : Value of the objective function at optimumCandidate.
%       optimizationData            : Struct containing detailed logs and process information.
%
%   Functionality:
%       - Computes gradients using finite differences.
%       - Searches along steepest descent direction.
%       - Handles design bounds by clipping.
%       - Modular and customizable architecture with pluggable functions.
%
%   Dependencies:
%       - evaluate_optimization_objective
%       - compute_gradient_finite_differences
%       - merge_name_value_pair_argument
%       - region_limit_line_search
%       - line search function (e.g., line_search_backtracking)
%       - convergence_criterion_optimum_candidate_variance (default)
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
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
    addParameter(p, 'EvaluationObjectiveOptions', {});
	addParameter(p, 'Alpha', 0.3, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Beta', 0.707, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
	parse(p, varargin{:});
	options = p.Results;

    [initialObjective,objectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction,design,options.EvaluationObjectiveOptions{:});

    stepSize = 1; % Initial step size
    stepSizeData = stepSize;  % Initialize step size data to store the values for each iteration
    stoppingCriterion = initialObjective - (stepSize/2) * (searchDirection * searchDirection');

    isOutputLinesearchData = (nargout>=3);
    if(isOutputLinesearchData)
        lineSearchData(1) = struct(...
            'StepSize',stepSize,...
            'CurrentObjective',initialObjective,...
            'ObjectiveEvaluationData',objectiveEvaluationData,...
            'StoppingCriterion',stoppingCriterion);
    end

    hasConverged = false;
    iteration = 1;
    while(~hasConverged && iteration < options.MaxIterations)
        nextDesign = design + stepSize * searchDirection;
        [nextObjective,nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});

        stoppingCriterion = initialObjective - (stepSize/2) * (searchDirection * searchDirection');
        hasConverged = (nextObjective <= stoppingCriterion);

        if(~hasConverged)
            stepSize = options.Beta * stepSize;
            stepSizeData = [stepSizeData, stepSize]; % Store the new step size for each iteration
        end

        if(isOutputLinesearchData)  
            lineSearchData(iteration+1) = struct(...
                'StepSize',stepSize,...
                'CurrentObjective',nextObjective,...
                'ObjectiveEvaluationData',nextObjectiveEvaluationData,...
                'StoppingCriterion',stoppingCriterion);
        end
        iteration = iteration + 1;
    end
end