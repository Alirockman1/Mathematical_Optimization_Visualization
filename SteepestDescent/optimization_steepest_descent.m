function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_steepest_descent(objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
%OPTIMIZATION_STEEPEST_DESCENT Perform unconstrained optimization using steepest descent with line search.
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE] = OPTIMIZATION_STEEPEST_DESCENT(OBJECTIVEFUNCTION, INITIALDESIGN,
%   DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND) performs optimization on the
%   given OBJECTIVEFUNCTION starting from INITIALDESIGN, constrained by the lower
%   and upper bounds of the design space, and returns the optimal solution and
%   its objective value.
%
%   [...] = OPTIMIZATION_STEEPEST_DESCENT(..., 'Name', Value, ...) accepts additional
%   name-value pair arguments to configure the optimization behavior:
%
%       'MaxIterations'              - Maximum number of iterations (default: 100)
%       'Tolerance'                 - Convergence tolerance on gradient norm (default: 1e-6)
%       'UseNormalizedGradientLineSearch' - Use unit-norm gradient for line search direction (default: false)
%       'LineSearchFunction'        - Function handle to custom line search method
%                                     (default: @line_search_backtracking)
%       'EvaluationObjectiveOptions'- Cell array of additional parameters for objective evaluation
%       'GradientOptions'           - Cell array of options passed to the gradient computation
%       'LineSearchOptions'         - Cell array of options passed to the line search function
%       'ConvergenceCriterionFunction' - Function handle for convergence check
%                                        (default: @convergence_criterion_optimum_candidate_variance)
%       'ConvergenceCriterionOptions' - Cell array of options for convergence criterion
%
%   OUTPUTS:
%       OPTIMUMCANDIDATE        - Optimal design point found by the algorithm
%       OPTIMUMOBJECTIVEVALUE   - Value of the objective function at OPTIMUMCANDIDATE
%       OPTIMIZATIONDATA        - Struct containing detailed data from each iteration:
%           .ProblemData
%           .InitialData
%           .FinalDesign
%           .FinalObjectiveValue
%           .FinalObjectiveValueData
%           .Iterations
%           .IterationData(iteration) - Struct array of iteration logs including:
%               - OptimumCandidate
%               - OptimumCandidateObjectiveValue
%               - Gradient and GradientData
%               - SearchDirection
%               - StepSize and StepSizeData
%               - Previous candidate and objective value
%
%   FUNCTIONALITY:
%       This implementation uses a basic steepest descent algorithm that iteratively
%       moves in the direction of the negative gradient of the objective function.
%       A custom line search function is used to determine the optimal step size
%       along the descent direction. The algorithm continues until convergence is
%       detected based on a user-defined criterion or the maximum number of iterations
%       is reached.
%
%   DEPENDENCIES:
%       - evaluate_optimization_objective
%       - compute_gradient_finite_differences
%       - region_limit_line_search
%       - line_search_backtracking (or custom)
%       - merge_name_value_pair_argument
%       - convergence_criterion_optimum_candidate_variance (or custom)
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
    addParameter(p, 'UseNormalizedGradientLineSearch', false, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'LineSearchFunction', @line_search_backtracking, @(x) isa(x, 'function_handle'));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'GradientOptions', {});
    addParameter(p, 'LineSearchOptions', {});
    addParameter(p,'ConvergenceCriterionFunction',@convergence_criterion_optimum_candidate_variance);
    addParameter(p,'ConvergenceCriterionOptions',{});
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
    [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction,currentDesign,options.EvaluationObjectiveOptions{:});

    optimizationData = struct();
    optimizationData.IterationData = repmat(... 
        struct('Iteration', [], ...
               'OptimumCandidate', [], ...
               'OptimumCandidateObjectiveValue', [], ...
               'OptimumCandidateObjectiveEvaluationData', [], ...
               'Gradient', [], ...
               'GradientData', [],...
               'StepSize', [], ...
               'StepSizeData', [], ...
               'SearchDirection', [], ...  
               'OptimumCandidatePrevious', [], ...  
               'OptimumCandidateObjectiveValuePrevious', []), ... 
        options.MaxIterations, 1);

    % only log if necessary
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
            'GradientData',[]);
    end

    hasConverged = false;

    iteration = 1;

    % Steepest descent loop
    while ~hasConverged && iteration <= options.MaxIterations
        % Compute gradient using numerical finite differences
        [gradient, gradientData] = compute_gradient_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, gradientOptions{:}); 
        
        % Compute steepest descent direction (normalized negative gradient)
        gradientNorm = norm(gradient);
        if gradientNorm < eps
            hasConverged = true;
        end

        searchDirection = -gradient;
        if(options.UseNormalizedGradientLineSearch)
            searchDirection = searchDirection / gradientNorm;
        end

        % Perform line search to determine step size
        maxStepSize = region_limit_line_search([],currentDesign,searchDirection,designSpace);
        [stepSize, stepSizeData] = options.LineSearchFunction(objectiveFunction, currentDesign, searchDirection, maxStepSize, lineSearchOptions{:});

        % Update design variables
        nextDesign = currentDesign + stepSize * searchDirection;
        nextDesign = max(min(nextDesign, designSpaceUpperBound), designSpaceLowerBound);

        % Evaluate objective function
        [nextObjectiveValue,nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction,nextDesign,options.EvaluationObjectiveOptions{:});

        % Convergence check
        hasConverged = options.ConvergenceCriterionFunction(...
            nextDesign,currentDesign,...
            nextObjectiveValue,currentObjectiveValue,...
            [],[],...
            [],[],...
            options.ConvergenceCriterionOptions{:});

        % Log iteration data
        optimizationData.IterationData(iteration) = struct('Iteration', iteration, ...
            'OptimumCandidate', nextDesign, ...
            'OptimumCandidateObjectiveValue', nextObjectiveValue, ...
            'OptimumCandidateObjectiveEvaluationData', nextObjectiveEvaluationData, ...
            'Gradient', gradient, ...
            'GradientData', gradientData, ...
            'StepSize', stepSize, ...
            'StepSizeData', stepSizeData, ...
            'SearchDirection', searchDirection, ...
            'OptimumCandidatePrevious', currentDesign, ...
            'OptimumCandidateObjectiveValuePrevious', currentObjectiveValue);

        % Update for next iteration
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

