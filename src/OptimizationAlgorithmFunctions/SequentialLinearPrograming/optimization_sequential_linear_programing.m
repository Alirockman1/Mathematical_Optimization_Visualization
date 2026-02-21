function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_sequential_linear_programing(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, moveLimit, designSpaceLowerBound, designSpaceUpperBound, varargin)
% FIND_BASIC_VARIABLES Identifies basic variables in a simplex tableau.
%
%   ISBASICVARIABLE = FIND_BASIC_VARIABLES(TABLEAU) determines which variables
%   in the simplex tableau are basic variables (variables in the basis).
%
%   INPUTS:
%       tableau - The simplex tableau matrix (m×n) where:
%                 - Rows 1:m-1 represent constraints
%                 - Row m represents the objective function
%                 - Columns 1:n-1 represent variables
%                 - Column n represents the RHS values
%
%   OUTPUTS:
%       isBasicVariable - Logical vector (1×n-1) where:
%                         - true indicates a basic variable
%                         - false indicates a non-basic variable
%
%   ALGORITHM:
%       1. Examines each variable column (excluding RHS column)
%       2. A variable is basic if its column contains:
%          a) Exactly one coefficient with absolute value > tolerance
%          b) The identified coefficient equals 1 within tolerance
%          c) All other coefficients in column are zero within tolerance
%       3. Uses tolerance of 1e-10 for numerical comparisons
%
%   IMPLEMENTATION NOTES:
%       - Handles numerical precision issues through tolerance
%       - Efficiently processes large tableaus
%       - Compatible with both standard and revised simplex methods
%
%   REFERENCES:
%       [1] Linear Programming and Network Flows, Bazaraa et al.
%       [2] Numerical Recipes: The Art of Scientific Computing, Press et al.
%
%   SEE ALSO:
%       solve_simplex_tableau, optimization_simplex, create_initial_tableau
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   SPDX-License-Identifier: Apache-2.0

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MinimumMoveLimit', 1e-4, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MaxFeasibilityRetries',10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SaveTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'PrintTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'GradientOptions', {});
    addParameter(p, 'ConvergenceCriterionFunction', @convergence_criterion_optimum_candidate_variance);
    addParameter(p, 'ConvergenceCriterionOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    currentDesign = initialDesign;
    currentMoveLimit = moveLimit;

    [initialObjectiveValue,initialObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, initialDesign, options.EvaluationObjectiveOptions{:});

    optimizationData = struct();
    optimizationData.IterationData = repmat(struct(...
            'Iteration', [], ...
            'PreviousOptimumCandidate', [], ...
            'PreviousOptimumCandidateObjectiveValue', [], ...
            'PreviousOptimumCandidateObjectiveEvaluationData', [], ...
            'PreviousMoveLimit', [], ...
            'OptimumCandidate', [], ...
            'OptimumCandidateObjectiveValue', [], ...
            'OptimumCandidateObjectiveEvaluationData', [], ...
            'MoveLimit', []), ...
        options.MaxIterations, 1);

    %DEL: only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
        optimizationData.ProblemData = struct(...
            'ObjectiveFunction',objectiveFunction,...
            'InequalityConstraintFunction',inequalityConstraintFunction,...
            'EqualityConstraintFunction',equalityConstraintFunction,...
            'DesignSpaceLowerBound',designSpaceLowerBound,...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options);

        optimizationData.InitialData = struct(...
            'OptimumCandidate',initialDesign,...
            'OptimumCandidateObjectiveValue',initialObjectiveValue,...
            'OptimumCandidateObjectiveEvaluationData',initialObjectiveEvaluationData);

    end

    iteration = 1;
    hasConverged = false;

    while ~hasConverged && iteration <= options.MaxIterations
        % Step 1: Linearize objective and constraints at current design point
        [gradientObjective, ~] = compute_gradient_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, options.GradientOptions{:});
        [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});

        objectiveFunctionLinearized = @(x) currentObjectiveValue + gradientObjective * (x(:) - currentDesign(:));

        if ~isempty(inequalityConstraintFunction)
            for i = 1:numel(inequalityConstraintFunction)
                [gradientConstraint, ~] = compute_gradient_finite_differences(inequalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.GradientOptions{:});
    
                constraintValue = inequalityConstraintFunction{i}(currentDesign);
    
                % Store the linearized constraint in the cell array
                inequalityConstraintFunctionLinearized{i} = @(x) constraintValue + gradientConstraint * (x(:) - currentDesign(:));
    
                % Extracting matrixes for test case
                A(i, :) = gradientConstraint(:)';
                b(i) = -1*(constraintValue + gradientConstraint(:)' * (-currentDesign(:)));
    
            end
        end

        if ~isempty(equalityConstraintFunction)
            for i = 1:numel(equalityConstraintFunction)
                [gradientConstraint, ~] = compute_gradient_finite_differences(equalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.GradientOptions{:});
    
                constraintValue = equalityConstraintFunction{i}(currentDesign);
    
                % Store the linearized constraint in the cell array
                equalityConstraintFunctionLinearized{i} = @(x) constraintValue + gradientConstraint * (x(:) - currentDesign(:));
    
                % Extracting matrixes for test case
                Aeq(i, :) = gradientConstraint(:)';
                beq(i) = -1*(constraintValue + gradientConstraint(:)' * (-currentDesign(:)));

            end
        else
            equalityConstraintFunctionLinearized = [];
        end

        % Extracting matrixes for test case
        f = gradientObjective(:)'; 

        lowerMoveLimit = currentDesign .* (1 - currentMoveLimit);
        upperMoveLimit = currentDesign .* (1 + currentMoveLimit);
        
        % Update bounds using element-wise min/max
        newDesignSpaceLowerBound = max(lowerMoveLimit, designSpaceLowerBound(:)');
        newDesignSpaceUpperBound = min(upperMoveLimit, designSpaceUpperBound(:)');

        % Step 2.0 === Compute the new design point ===
       [nextDesign, ~, simplexOptimizationData] = optimization_simplex(objectiveFunctionLinearized, inequalityConstraintFunction, equalityConstraintFunction, currentDesign,...
           newDesignSpaceLowerBound, newDesignSpaceUpperBound,'MaxIterations', 30, 'Tolerance', 1e-6);
       nextDesign = nextDesign';
   
        nextMoveLimit = max(0.9 * currentMoveLimit, options.MinimumMoveLimit);

        % Step 2.1 === Check Feasibility ===
        feasible = true;
        tolerance = options.Tolerance;

        % Check inequality constraints
        for i = 1:length(inequalityConstraintFunction)
            constraintVal = inequalityConstraintFunction{i}(nextDesign);
            if any(constraintVal > tolerance)
                feasible = false;
                break;
            end
        end

        % Check bounds
        if any(nextDesign < designSpaceLowerBound - tolerance) || any(nextDesign > designSpaceUpperBound + tolerance)
            feasible = false;
        end

        % Retry with reduced move limits until a feasible design is found
        maxFeasibilityRetries = options.MaxFeasibilityRetries;
        feasibilityRetryCount = 0;

        
        while ~feasible && feasibilityRetryCount < maxFeasibilityRetries
            fprintf('Iteration %d: Infeasible design. Reducing move limit to %.4f and retrying...\n', iteration, nextMoveLimit);
        
            % Feasible design
            feasibleDesign = currentDesign;
            feasibleMoveLimit = nextMoveLimit;
            
            % Recompute move limits and bounds
            lowerMoveLimit = feasibleDesign .* (1 - feasibleMoveLimit);
            upperMoveLimit = feasibleDesign .* (1 + feasibleMoveLimit);
            newDesignSpaceLowerBound = max(lowerMoveLimit, designSpaceLowerBound(:)');
            newDesignSpaceUpperBound = min(upperMoveLimit, designSpaceUpperBound(:)');
        
            % Recalculate deltaDesign with reduced move limits
            [nextDesign, ~, ~] = optimization_simplex(objectiveFunctionLinearized, inequalityConstraintFunction, equalityConstraintFunction, currentDesign,...
                newDesignSpaceLowerBound, newDesignSpaceUpperBound,'MaxIterations', 30, 'Tolerance', 1e-6);
            nextDesign = nextDesign';

            nextMoveLimit = max(0.9 * feasibleMoveLimit, options.MinimumMoveLimit);
        
            % Constraint check
            for i = 1:length(inequalityConstraintFunction)
                constraintVal = inequalityConstraintFunction{i}(nextDesign);
                if any(constraintVal > tolerance)
                    feasible = false;
                    break;
                end
            end
        
            % Bounds check
            if any(nextDesign < designSpaceLowerBound - tolerance) || any(nextDesign > designSpaceUpperBound + tolerance)
                feasible = false;
            end
        
            feasibilityRetryCount = feasibilityRetryCount + 1;
        end

        [nextObjectiveValue,nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});


        % Step 3 === Convergence Check ===
        hasConverged = options.ConvergenceCriterionFunction(...
        nextDesign, currentDesign, ...
        nextObjectiveValue, currentObjectiveValue, ...
        [], [], [], [], ...
        options.ConvergenceCriterionOptions{:});

        % Log iteration data
        optimizationData.IterationData(iteration) = struct(...
            'Iteration', iteration, ...
            'PreviousOptimumCandidate', currentDesign, ...
            'PreviousOptimumCandidateObjectiveValue', currentObjectiveValue, ...
            'PreviousOptimumCandidateObjectiveEvaluationData', currentObjectiveEvaluationData, ...
            'PreviousMoveLimit', currentMoveLimit, ...
            'OptimumCandidate', nextDesign, ...
            'OptimumCandidateObjectiveValue', nextObjectiveValue, ...
            'OptimumCandidateObjectiveEvaluationData', nextObjectiveEvaluationData, ...
            'MoveLimit', nextMoveLimit);

        iteration = iteration +1;

        % Update the move limit and design point
        currentMoveLimit = nextMoveLimit;
        currentDesign    = nextDesign;

    end

    % Trim unused iterations from the struct array
    optimizationData.IterationData = optimizationData.IterationData(1:iteration - 1);

    % Outputs
    optimumCandidate = nextDesign;
    optimumObjectiveValue = nextObjectiveValue;

    % Final log
    optimizationData.FinalDesign = optimumCandidate;
    optimizationData.FinalObjectiveValue = optimumObjectiveValue;
    optimizationData.FinalObjectiveValueData = nextObjectiveEvaluationData;
    optimizationData.Iterations = iteration - 1;
end
