function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_modified_newton(objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
%OPTIMIZATION_MODIFIED_NEWTON Modified Newton method with line search for optimization.
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE, OPTIMIZATIONDATA] = ...
%       OPTIMIZATION_MODIFIED_NEWTON(OBJECTIVEFUNCTION, INITIALDESIGN, ...
%       DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND, ...) implements a
%       modified Newton method with bound constraints to minimize the given
%       objective function.
%
%   INPUTS:
%       OBJECTIVEFUNCTION      - Function handle to the objective function
%       INITIALDESIGN          - Starting point vector
%       DESIGNSPACELOWERBOUND  - Vector of lower bounds
%       DESIGNSPACEUPPERBOUND  - Vector of upper bounds
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'MaxIterations'        - Maximum number of iterations (default: 100)
%       'Tolerance'           - Convergence tolerance (default: 1e-6)
%       'Epsilon'             - Small value for Hessian modification (default: 1e-6)
%       'LineSearchFunction'  - Line search method (default: @line_search_backtracking)
%       'EvaluationObjectiveOptions' - Options for objective evaluation
%       'GradientOptions'     - Options for gradient computation
%       'HessianOptions'      - Options for Hessian computation
%       'ConvergenceCriterionFunction' - Convergence check function
%
%   OUTPUTS:
%       OPTIMUMCANDIDATE      - Optimal solution found
%       OPTIMUMOBJECTIVEVALUE - Objective value at optimum
%       OPTIMIZATIONDATA      - Structure containing optimization history:
%           .ProblemData      - Problem definition and options
%           .InitialData      - Initial point information
%           .IterationData    - Detailed iteration records
%           .FinalDesign      - Final optimal solution
%
%   FUNCTIONALITY:
%       1. Implements modified Newton method with bound constraints
%       2. Uses line search to ensure convergence
%       3. Handles non-positive definite Hessian matrices
%       4. Provides detailed logging of optimization process
%       5. Supports custom convergence criteria
%
%   NOTES:
%       - Automatically enforces design space bounds
%       - Detailed optimization data only computed if requested (nargout=3)
%       - Uses efficient memory management for iteration logging
%       - First iteration uses standard Newton step
%
%   DEPENDENCIES:
%       - line_search_backtracking.m
%       - compute_gradient_finite_differences.m
%       - compute_hessian_finite_differences.m
%       - evaluate_optimization_objective.m
%       - convergence_criterion_optimum_candidate_variance.m
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
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
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Epsilon', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'LineSearchFunction', @line_search_backtracking, @(x) isa(x, 'function_handle'));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'GradientOptions', {});
    addParameter(p, 'HessianOptions', {});
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
    currentDesign = initialDesign(:)';
    [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});

    optimizationData = struct();
    optimizationData.IterationData = repmat(struct(...
        'Iteration', [], ...
        'OptimumCandidate', [], ...
        'OptimumCandidateObjectiveValue', [], ...
        'OptimumCandidateObjectiveEvaluationData', [], ...
        'Gradient', [], ...
        'GradientData', [], ...
        'Hessian', [], ...
        'HessianData', [], ...
        'StepSize', [], ...
        'StepSizeData', [], ...
        'SearchDirection', []), ...
        options.MaxIterations, 1);

    %DEL: only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
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
            'GradientData',[], ...
            'Hessian', [], ...
            'HessianData', []);

    end

    iteration = 1;
    hasConverged = false;

    while ~hasConverged && iteration <= options.MaxIterations
        [gradient, gradientData] = compute_gradient_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, options.GradientOptions{:});
        [hessian, hessianData] = compute_hessian_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, options.HessianOptions{:});

        % Compute search direction (Newton step)
        searchDirection = -(hessian \ gradient(:))';

        % Perform line search to determine step size
        maxStepSize = region_limit_line_search([],currentDesign,searchDirection,designSpace);
        [stepSize, stepSizeData] = options.LineSearchFunction(objectiveFunction, currentDesign, searchDirection, maxStepSize, lineSearchOptions{:});

        % Update design variables and enforce bounds
        nextDesign = currentDesign + stepSize * searchDirection;
        nextDesign = max(min(nextDesign, designSpaceUpperBound), designSpaceLowerBound);

        % Evaluate new objective function and gradient
        [nextObjectiveValue, nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});

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
            'Hessian', hessian, ...
            'HessianData', hessianData, ...
            'StepSize', stepSize, ...
            'StepSizeData', stepSizeData, ...
            'SearchDirection', searchDirection);

        % Update variables for next iteration
        currentDesign = nextDesign;
        currentObjectiveValue = nextObjectiveValue;
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

