function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_conjugate_gradient(objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
%OPTIMIZATION_CONJUGATE_GRADIENT Conjugate gradient optimization with line search.
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE, OPTIMIZATIONDATA] = ...
%       OPTIMIZATION_CONJUGATE_GRADIENT(OBJECTIVEFUNCTION, INITIALDESIGN, ...
%       DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND, ...) implements the
%       conjugate gradient method to minimize OBJECTIVEFUNCTION within specified
%       bounds, returning the optimal solution and optimization history.
%
%   INPUTS:
%       OBJECTIVEFUNCTION      - Function handle to objective (must return scalar)
%       INITIALDESIGN          - Starting point (column vector)
%       DESIGNSPACELOWERBOUND  - Vector of lower bounds
%       DESIGNSPACEUPPERBOUND  - Vector of upper bounds
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'MaxIterations'        - Maximum iterations (default: 100)
%       'Tolerance'           - Convergence threshold (default: 1e-6)
%       'LineSearchFunction'  - Line search method (default: @line_search_backtracking)
%       'BetaFunction'        - Beta computation method (default: @compute_beta_fletcher_reeves)
%       'EvaluationObjectiveOptions' - Options for objective evaluation
%       'GradientOptions'     - Options for gradient computation
%       'LineSearchOptions'   - Options for line search
%       'ConvergenceCriterionFunction' - Convergence check (default: @convergence_criterion_optimum_candidate_variance)
%
%   OUTPUTS:
%       OPTIMUMCANDIDATE      - Optimal solution found
%       OPTIMUMOBJECTIVEVALUE - Objective value at optimum
%       OPTIMIZATIONDATA      - Structure containing complete optimization history:
%           .ProblemData      - Problem definition and options
%           .InitialData      - Initial point information
%           .IterationData    - Array of iteration records
%           .FinalDesign      - Final optimal solution
%
%   FUNCTIONALITY:
%       1. Implements nonlinear conjugate gradient method
%       2. Supports various beta computation methods (Fletcher-Reeves default)
%       3. Performs bound-constrained optimization
%       4. Provides detailed iteration logging
%       5. Supports custom convergence criteria
%
%   NOTES:
%       - First iteration uses steepest descent
%       - Automatically handles bound constraints during line search
%       - Detailed optimization data only computed if requested (nargout=3)
%       - Uses efficient memory management for iteration logging
%
%   DEPENDENCIES:
%       - line_search_backtracking.m
%       - compute_gradient_finite_differences.m
%       - compute_beta_fletcher_reeves.m
%       - convergence_criterion_optimum_candidate_variance.m
%       - evaluate_optimization_objective.m
%       - merge_name_value_pair_argument.m
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
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

   % Parse input arguments
    % Parse input arguments
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'LineSearchFunction', @line_search_backtracking, @(x) isa(x, 'function_handle'));
    addParameter(p, 'BetaFunction', @compute_beta_fletcher_reeves, @(x) isa(x, 'function_handle'));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'GradientOptions', {});
    addParameter(p, 'LineSearchOptions', {});
    addParameter(p, 'ConvergenceCriterionFunction', @convergence_criterion_optimum_candidate_variance);
    addParameter(p, 'ConvergenceCriterionOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    % Validate inputs
    if length(initialDesign) ~= length(designSpaceLowerBound) || ...
       length(initialDesign) ~= length(designSpaceUpperBound)
        error('Initial design and design space bounds must have the same dimensions.');
    end

    % merge options
    [~,gradientOptions] = merge_name_value_pair_argument(options.GradientOptions,{'EvaluationObjectiveOptions',options.EvaluationObjectiveOptions});
    [~,lineSearchOptions] = merge_name_value_pair_argument(options.LineSearchOptions,{'EvaluationObjectiveOptions',options.EvaluationObjectiveOptions});
    designSpace = [designSpaceLowerBound;designSpaceUpperBound];

    % Initialization
    currentDesign = initialDesign(:)';  % Ensure row vector
    [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});

    % only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
        optimizationData = struct();

        optimizationData.ProblemData = struct(...
            'ObjectiveFunction',objectiveFunction,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options);

        optimizationData.InitialData = struct(...
            'OptimumCandidate',currentDesign,...
            'OptimumCandidateObjectiveValue',currentObjectiveValue,...
            'OptimumCandidateObjectiveEvaluationData',currentObjectiveEvaluationData,...
            'Gradient',[],...
            'GradientData',[]);

            optimizationData.IterationData = repmat(struct(...
                'Iteration', [], ...
                'OptimumCandidate', [], ...
                'OptimumCandidateObjectiveValue', [], ...
                'OptimumCandidateObjectiveEvaluationData', [], ...
                'Gradient', [], ...
                'GradientData', [], ...
                'StepSize', [], ...
                'StepSizeData', [], ...
                'SearchDirection', [], ...
                'Beta', []), ...
                options.MaxIterations, 1);
    end

    % First iteration uses Steepest Descent
    gradientPrevious = zeros(1,size(currentDesign,2));
    searchDirectionPrevious = zeros(1,size(currentDesign,2));

    hasConverged = false;
    iteration = 1;

    % Conjugate Gradient loop
    while ~hasConverged && iteration <= options.MaxIterations
        % Compute gradient
        [gradient, gradientData] = compute_gradient_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, gradientOptions{:});

        % Compute Fletcher-Reeves beta
        beta = options.BetaFunction(gradient, gradientPrevious);

        % Compute new search direction
        searchDirection = -gradient + beta * searchDirectionPrevious;

        % Perform line search to determine step size
        maxStepSize = region_limit_line_search([],currentDesign,searchDirection,designSpace);
        [stepSize, stepSizeData] = options.LineSearchFunction(objectiveFunction, currentDesign, searchDirection, maxStepSize, lineSearchOptions{:});

        % Update design variables
        nextDesign = currentDesign + stepSize * searchDirection;
        nextDesign = max(min(nextDesign, designSpaceUpperBound), designSpaceLowerBound);

        % Evaluate objective function
        [nextObjectiveValue,nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});

        % Convergence check
        hasConverged = options.ConvergenceCriterionFunction(...
            nextDesign, currentDesign, ...
            nextObjectiveValue, currentObjectiveValue, ...
            [], [], [], [], ...
            options.ConvergenceCriterionOptions{:});

        % Log iteration data
        optimizationData.IterationData(iteration) = struct(...
            'Iteration', iteration, ...
            'OptimumCandidate', nextDesign, ...
            'OptimumCandidateObjectiveValue', nextObjectiveValue, ...
            'OptimumCandidateObjectiveEvaluationData', nextObjectiveEvaluationData, ...
            'Gradient', gradient, ...
            'GradientData', gradientData, ...
            'StepSize', stepSize, ...
            'StepSizeData', stepSizeData, ...
            'SearchDirection', searchDirection, ...
            'Beta', beta);

        % Update for next iteration
        currentDesign = nextDesign;
        currentObjectiveValue = nextObjectiveValue;
        gradientPrevious = gradient;
        searchDirectionPrevious = searchDirection;
        iteration = iteration + 1;
    end

    % Trim unused iterations from the struct array
    optimizationData.IterationData = optimizationData.IterationData(1:iteration - 1);

    % Outputs
    optimumCandidate = currentDesign;
    optimumObjectiveValue = currentObjectiveValue;

    % Final log
    optimizationData.FinalDesign = optimumCandidate;
    optimizationData.FinalObjectiveValue = optimumObjectiveValue;
    optimizationData.FinalObjectiveValueData = currentObjectiveEvaluationData;
    optimizationData.Iterations = iteration - 1;
end

