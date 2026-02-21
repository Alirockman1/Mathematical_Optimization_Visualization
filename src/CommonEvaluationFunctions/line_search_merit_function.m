function [optimalStepSize, stepSizeData, linesearchData] = line_search_merit_function(objectiveFunction, equalityConstraintFunctions, inequalityConstraintFunctions, ...
    currentDesign, searchDirection, previousWeights, varargin)
%LINE_SEARCH_MERIT_FUNCTION Merit function-based line search for constrained optimization.
%
%   [OPTIMALSTEPSIZE, STEPSIZEDATA, LINESEARCHDATA] = LINE_SEARCH_MERIT_FUNCTION(OBJECTIVEFUNCTION,
%   EQUALITYCONSTRAINTFUNCTIONS, INEQUALITYCONSTRAINTFUNCTIONS, CURRENTDESIGN,
%   SEARCHDIRECTION, PREVIOUSWEIGHTS, ...) performs a merit function-based line search
%   for constrained optimization problems by evaluating merit function values across
%   a range of step sizes to find the optimal step.
%
%   INPUTS:
%       OBJECTIVEFUNCTION              - Function handle to evaluate objective
%       EQUALITYCONSTRAINTFUNCTIONS    - Cell array of equality constraint functions
%       INEQUALITYCONSTRAINTFUNCTIONS  - Cell array of inequality constraint functions
%       CURRENTDESIGN                  - Current design point
%       SEARCHDIRECTION                - Direction for line search
%       PREVIOUSWEIGHTS                - Structure with penalty weights:
%                                       .equality - weights for equality constraints
%                                       .inequality - weights for inequality constraints
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'Alpha'                        - Convergence parameter (default: 0.8)
%       'Beta'                         - Backtracking parameter (default: 0.707)
%       'MaxIterations'                - Maximum iterations (default: 100)
%       'EvaluationObjectiveOptions'   - Options for objective evaluation (default: {})
%
%   OUTPUTS:
%       OPTIMALSTEPSIZE               - Optimal step size found
%       STEPSIZEDATA                  - Structure containing detailed search data:
%           .StepSize                 - Array of tested step sizes
%           .DesignVariable           - Matrix of design points tested
%           .ObjectiveValue           - Array of objective values
%           .MeritValue               - Array of merit function values
%           .OptimalStepSize          - Array of optimal step sizes found
%           .MinimumMerit             - Array of minimum merit values
%       LINESEARCHDATA                - Structure containing iteration history:
%           .InitialData              - Initial evaluation data
%           .IterationData            - Array of iteration-specific data
%
%   FUNCTIONALITY:
%       1. Evaluates merit function across uniform step size range [0, 1]
%       2. Merit function combines objective and constraint violations
%       3. Finds step size that minimizes merit function value
%       4. Uses weighted penalty approach for constraint handling
%       5. Provides convergence check based on merit improvement
%       6. Maintains complete search history for analysis
%
%   MERIT FUNCTION FORMULATION:
%       Merit(x) = f(x) + sum(w_eq * |h_i(x)|) + sum(w_ineq * max(0, -g_j(x)))
%       where:
%           f(x)    - objective function
%           h_i(x)  - equality constraints
%           g_j(x)  - inequality constraints (g_j(x) >= 0)
%           w_eq    - equality constraint weights
%           w_ineq  - inequality constraint weights
%
%   NOTES:
%       - Uses uniform sampling of step sizes from 0 to 1
%       - Merit function handles both equality and inequality constraints
%       - Convergence based on sufficient decrease condition
%       - Detailed logging available for analysis and debugging
%       - Handles empty constraint sets gracefully
%
%   REFERENCES:
%       [1] Nocedal, J. & Wright, S.J. (2006). "Numerical Optimization", 2nd ed.
%       [2] Fletcher, R. (1987). "Practical Methods of Optimization", 2nd ed.
%
%   SEE ALSO:
%       line_search_golden_ratio, line_search_backtracking, 
%       evaluate_optimization_objective, optimization_barrier_penalty
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

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'Alpha', 0.8, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Beta', 0.707, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'MaxIterations', 100, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    optimizationData = struct();
    linesearchData.IterationData = repmat(struct(...
        'Iteration', [], ...
        'StepSize', [],...
        'MeritValue', [], ...
        'Objective', [],...
        'EqualityConstraints', [],...
        'InequalityConstraints', [],...
        'StepSizeData', []), ...
        options.MaxIterations, 1);

    stepSize = 1;
    stepSizeData = stepSize;

    % Evaluate objective at currentDesign
    [initialObjective, ~] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});

    % Evaluate constraints at currentDesign
    initialEquality = evaluate_constraints(equalityConstraintFunctions, currentDesign);
    initialInequality = evaluate_constraints(inequalityConstraintFunctions, currentDesign);

    % Compute initial merit
    initialMerit = compute_merit(initialObjective, initialEquality, initialInequality, previousWeights);

    %DEL: only log if necessary
    isOutputLinesearchData = (nargout>=3);
    if(isOutputLinesearchData)
        linesearchData.InitialData = struct(...
            'StepSize', stepSize,...
            'MeritValue', initialMerit, ...
            'Objective', initialObjective,...
            'EqualityConstraints', initialEquality,...
            'InequalityConstraints', initialInequality);
    end

    iterationCount = 1;
    hasConverged = false;

    while ~hasConverged && iterationCount < options.MaxIterations
        % Define a search range for stepSize (alpha) between 0 and 1
        stepSizeRange = linspace(0, 1, 100);  % 100 values between 0 and 1 for alpha (step size)
        minimumMerit = Inf;  % Initialize the minimum merit to a large value
        stepSizeData = struct('StepSize', [],...
            'DesignVariable', [],...
            'ObjectiveValue', [],...
            'MeritValue', [],...
            'OptimalStepSize',[],...
            'MinimumMerit',[]);
    
        % Try different values of alpha (step size) and find the one that minimizes the merit function
        for i = 1:length(stepSizeRange)
            % Update the design using the candidate step size (alpha)
            nextDesign = currentDesign + stepSizeRange(i) * searchDirection;
            nextObjective = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});
            nextEquality = evaluate_constraints(equalityConstraintFunctions, nextDesign);
            nextInequality = evaluate_constraints(inequalityConstraintFunctions, nextDesign);
    
            % Compute the merit value at the new design
            nextMerit = compute_merit(nextObjective, nextEquality, nextInequality, previousWeights);
    
            % Store the current alpha (step size), design, objective, and merit
            stepSizeData.StepSize = [stepSizeData.StepSize, stepSizeRange(i)];
            stepSizeData.DesignVariable = [stepSizeData.DesignVariable; nextDesign];
            stepSizeData.ObjectiveValue = [stepSizeData.ObjectiveValue, nextObjective];
            stepSizeData.MeritValue = [stepSizeData.MeritValue, nextMerit];
    
            % Update the minimum merit and optimal alpha (step size)
            if nextMerit < minimumMerit
                minimumMerit = nextMerit;
                optimalStepSize = stepSizeRange(i);
                stepSizeData.OptimalStepSize = [stepSizeData.OptimalStepSize, optimalStepSize];
                stepSizeData.MinimumMerit = [stepSizeData.MinimumMerit, minimumMerit];
            end
        end
    
        % Update the main step size with the optimal alpha (step size)
        optimalDesign = currentDesign + optimalStepSize * searchDirection;
        optimalObjectiveValue = evaluate_optimization_objective(objectiveFunction, optimalDesign, options.EvaluationObjectiveOptions{:});
        optimalEqualityValue = evaluate_constraints(equalityConstraintFunctions, optimalDesign);
        optimalInequalityValue = evaluate_constraints(inequalityConstraintFunctions, optimalDesign);
        optimalMeritValue = compute_merit(optimalObjectiveValue, optimalEqualityValue, optimalInequalityValue, previousWeights);
    
        % Check for convergence
        hasConverged = optimalMeritValue <= initialMerit - options.Alpha * stepSize * (searchDirection * searchDirection');
    
        % Store the final step size, merit, and relevant data after the loop
        if isOutputLinesearchData
            linesearchData.IterationData = struct('StepSize', stepSize,...
                'MeritValue', optimalMeritValue, ...
                'Objective', optimalObjectiveValue,...
                'EqualityConstraints', optimalEqualityValue,...
                'InequalityConstraints', optimalInequalityValue,...
                'StepSizeData', stepSizeData);  % Include the step size data
        end

        % Update iteration count
        iterationCount = iterationCount + 1;
    end
    
    % Return the final step size and the step size data
    optimalStepSize = stepSize;  % This is the optimal step size found in the loop
    % stepSizeData.OptimalAlpha = optimalAlpha;  % Store the optimal alpha (step size) found
    % stepSizeData.FinalMerit = finalMerit;  % Store the final merit value

% Return these values to the main function
end

    function constraintValues = evaluate_constraints(constraintFunctions, designPoint)
        if isempty(constraintFunctions)
            constraintValues = [];
        else
            constraintValues = zeros(length(constraintFunctions), 1);
            for idx = 1:length(constraintFunctions)
                constraintValues(idx) = constraintFunctions{idx}(designPoint);
            end
        end
    end

    function phi = compute_merit(f, h, g, weights)
        phi = f;
        if ~isempty(h)
            phi = phi + sum(weights.equality .* abs(h));
        end
        if ~isempty(g)
            phi = phi + sum(weights.inequality .* max(0, -g));
        end
    end