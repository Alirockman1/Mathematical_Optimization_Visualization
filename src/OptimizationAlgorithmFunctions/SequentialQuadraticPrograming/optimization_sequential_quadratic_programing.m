function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_sequential_quadratic_programing(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
% OPTIMIZATION_SEQUENTIAL_QUADRATIC_PROGRAMMING Solves nonlinear optimization problems using Sequential Quadratic Programming (SQP)
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE] = OPTIMIZATION_SEQUENTIAL_QUADRATIC_PROGRAMMING(OBJECTIVEFUNCTION, INEQUALITYCONSTRAINTFUNCTION, 
%   EQUALITYCONSTRAINTFUNCTION, INITIALDESIGN, DESIGNSPACELOWERBOUND, DESIGNSPACEUPPERBOUND) solves constrained nonlinear optimization problems 
%   by iteratively solving quadratic subproblems.
%
%   INPUTS:
%       objectiveFunction - Function handle for the objective function to minimize
%       inequalityConstraintFunction - Cell array of function handles for inequality constraints (g(x) ≤ 0)
%       equalityConstraintFunction - Cell array of function handles for equality constraints (h(x) = 0)
%       initialDesign - Initial guess for the design variables
%       designSpaceLowerBound - Vector of lower bounds for design variables
%       designSpaceUpperBound - Vector of upper bounds for design variables
%
%   OPTIONAL PARAMETERS (name-value pairs):
%       'MaxIterations' - Maximum number of iterations (default: 100)
%       'Tolerance' - Convergence tolerance (default: 1e-10)
%       'LagrangeTolerance' - Tolerance for Lagrange multiplier convergence (default: 1e-10)
%       'LineSearchFunction' - Function handle for line search method (default: @line_search_merit_function)
%       'InitialEqualityLagrangeMultiplier' - Initial values for equality constraint multipliers
%       'InitialInequalityLagrangeMultiplier' - Initial values for inequality constraint multipliers
%       'EvaluationObjectiveOptions' - Additional options for objective function evaluation
%       'GradientOptions' - Options for gradient computation
%       'HessianOptions' - Options for Hessian computation
%       'LineSearchOptions' - Options for line search
%       'ConvergenceCriterionFunction' - Function handle for convergence checking
%       'ConvergenceCriterionOptions' - Options for convergence criterion
%
%   OUTPUTS:
%       optimumCandidate - Optimal design point found
%       optimumObjectiveValue - Objective function value at optimum
%       optimizationData - Structure containing detailed optimization data
%
%   ALGORITHM:
%       Implements Sequential Quadratic Programming (SQP) method:
%       1. At each iteration, constructs a quadratic approximation of the Lagrangian
%       2. Solves the quadratic subproblem using active-set method
%       3. Performs line search using merit function
%       4. Updates design variables and Lagrange multipliers
%       5. Repeats until convergence criteria are met
%
%   IMPLEMENTATION DETAILS:
%       - Uses finite differences for gradient/Hessian if analytical forms not provided
%       - Handles both equality and inequality constraints
%       - Implements bound constraints on design variables
%       - Includes comprehensive logging of optimization progress
%
%   REFERENCES:
%       [1] Nocedal, J., & Wright, S. J. (2006). Numerical Optimization.
%       [2] Boggs, P. T., & Tolle, J. W. (1995). Sequential Quadratic Programming.
%
%   SEE ALSO:
%       optimization_active_set, line_search_merit_function, optimization_simplex
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   SPDX-License-Identifier: Apache-2.0

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'LagrangeTolerance', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'LineSearchFunction', @line_search_merit_function, @(x) isa(x, 'function_handle'));
    addParameter(p, 'InitialEqualityLagrangeMultiplier', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'InitialInequalityLagrangeMultiplier', [], @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'GradientOptions', {});
    addParameter(p, 'HessianOptions', {});
    addParameter(p, 'LineSearchOptions', {});
    addParameter(p, 'ConvergenceCriterionFunction', @convergence_criterion_optimum_candidate_variance);
    addParameter(p, 'ConvergenceCriterionOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    % Number of dimensions
    numberDimensions = length(initialDesign);

    if ~isempty(equalityConstraintFunction)
        numberEqualityConstraints = length(equalityConstraintFunction(initialDesign));
    else
        numberEqualityConstraints = 0;
    end

    % Initial lagrange multipliers
    if isempty(options.InitialEqualityLagrangeMultiplier)
        initialEqualityLagrangeMultiplier = zeros(numel(equalityConstraintFunction),1);
    else
        initialEqualityLagrangeMultiplier = options.InitialEqualityLagrangeMultiplier;
    end

    if isempty(options.InitialInequalityLagrangeMultiplier)
        initialInequalityLagrangeMultiplier = zeros(numel(inequalityConstraintFunction),1);
    else
        initialInequalityLagrangeMultiplier = options.InitialInequalityLagrangeMultiplier;
    end

    % Evaluate the objective function
    [initialObjectiveValue,initialObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, initialDesign, options.EvaluationObjectiveOptions{:});

    nextObjectiveEvaluationData = struct( ...
        'DesignPoint', [], ...
        'ObjectiveValueBase', [], ...
        'ObjectiveValueUnconstrained', []);

    stepSizeData = struct( ...
        'StepSize', [], ...
        'DesignVariable', [], ...
        'ObectiveValue', [], ...
        'MeritValue', [], ...
        'OptimalStepSize', [], ...
        'MinimumMerit', []);

    optimizationData = struct();
    optimizationData.IterationData = repmat(struct(...
        'Iteration', [], ...
        'OptimumCandidate', [], ...
        'OptimumCandidateObjectiveValue', [], ...
        'OptimumCandidateObjectiveEvaluationData', nextObjectiveEvaluationData, ...
        'GradientObjectiveFunction', [], ...
        'HessianObjectiveFunction', [], ...
        'EqualityLagrangeMultiplier', [], ...
        'InequalityLagrangeMultiplier', [], ...
        'MetaFunctionEqualityWeights', [],...
        'MetaFunctionInequalityWeights', [],...
        'StepSize', [], ...
        'StepSizeData', stepSizeData, ...
        'SearchDirection', [], ...
        'OptimumCandidatePrevious', [], ...
        'OptimumCandidateObjectiveValuePrevious', []), ...
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

    % Initializing parameters for the optimization loop
    iteration = 1;
    hasConverged = false;
    currentDesign = initialDesign;
    currentEqualityLagrangeMultiplier = initialEqualityLagrangeMultiplier;
    currentInequalityLagrangeMultiplier = initialInequalityLagrangeMultiplier;
    currentWeights.equality = zeros(numel(equalityConstraintFunction), 1);

    while ~hasConverged && iteration <= options.MaxIterations

        % Step 1: Derive Quadratic Sub Problem

        % Step 1-a: Compute gradients and values of objective function
        [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});
        [gradientObjective, ~] = compute_gradient_finite_differences(objectiveFunction, currentDesign,  'IsObjectiveFunction', true, options.GradientOptions{:});
        [hessianObjective, ~]  = compute_hessian_finite_differences(objectiveFunction, currentDesign, 'IsObjectiveFunction', true, options.HessianOptions{:});

        % Step 1-b: Compute gradients and values of constraints (if any)
        if ~isempty(equalityConstraintFunction)
            gradientEqualityConstraint = zeros(numel(equalityConstraintFunction),numberDimensions);
            hessianEqualityConstraint = zeros(numberDimensions);

            for i = 1:numel(equalityConstraintFunction)
                [gradientConstraint, ~] = compute_gradient_finite_differences(equalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.GradientOptions{:});
                [hessianConstraint, ~] = compute_hessian_finite_differences(equalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.HessianOptions{:});
    
                valueEqualityConstraint(i,1) = equalityConstraintFunction{i}(currentDesign);

                gradientEqualityConstraint(i,:) = gradientConstraint;
                hessianEqualityConstraint = hessianEqualityConstraint + currentEqualityLagrangeMultiplier(i).*hessianConstraint;

            end
        else
            valueEqualityConstraint    = [];
        end
    
        if ~isempty(inequalityConstraintFunction)
            gradientInequalityConstraint = zeros(numel(inequalityConstraintFunction),numberDimensions);
            hessianInequalityConstraint = zeros(numberDimensions);

            for i = 1:numel(inequalityConstraintFunction)
                [gradientConstraint, ~] = compute_gradient_finite_differences(inequalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.GradientOptions{:});
                [hessianConstraint, ~] = compute_hessian_finite_differences(inequalityConstraintFunction{i}, currentDesign, 'IsObjectiveFunction', false, options.HessianOptions{:});

                valueInequalityConstraint(i,1) = inequalityConstraintFunction{i}(currentDesign);

                gradientInequalityConstraint(i,:) = gradientConstraint;
                hessianInequalityConstraint =  hessianInequalityConstraint +  currentInequalityLagrangeMultiplier(i).*hessianConstraint;

            end
        else
            valueInequalityConstraint    = [];
        end
    
        % Step 1-c: Create the quadratic problem's matrices
        % Creating the Q Matrix
        if ~isempty(inequalityConstraintFunction) && ~isempty(equalityConstraintFunction)
            objective.Q = hessianObjective + hessianInequalityConstraint + hessianEqualityConstraint;
        elseif ~isempty(inequalityConstraintFunction)
            objective.Q = hessianObjective + hessianInequalityConstraint;
        elseif ~isempty(equalityConstraintFunction)
            objective.Q = hessianObjective + hessianEqualityConstraint;
        else
            objective.Q = hessianObjective;
        end
    
        % Creating the e vector
        objective.e = gradientObjective;
    
        % Creating the A matrix and b vector
        if ~isempty(equalityConstraintFunction)
            A = gradientEqualityConstraint;
            b = -valueEqualityConstrain;
            equalityConstraint.A = A;                                      % Gradient of the equality constraints
            equalityConstraint.b = b;                                      % Values of the equality constraints
        else
            equalityConstraint.A = [];
            equalityConstraint.b = [];
        end
        
        % Creating the C matrix and d vector
        if ~isempty(inequalityConstraintFunction)
            C = gradientInequalityConstraint;
            d = -valueInequalityConstraint;
            inequalityConstraint.C = C;                                    % Gradient of the equality constraints
            inequalityConstraint.d = d;                                    % Values of the equality constraints
        else
            inequalityConstraint.C = [];
            inequalityConstraint.d = [];
        end
    
        % Step 2: Solve the quadratic sub-problem using quadprog (Linear Programming solver)
        % Formulate the problem for linprog: min (1/2 * x'Qx + e'x) subject to Ax = b and A_ineqx <= b_ineq

        [activeSetDesign, ~, activeSetData] = optimization_active_set(objectiveFunction, inequalityConstraintFunction,...
            equalityConstraintFunction, currentDesign, designSpaceLowerBound, designSpaceUpperBound);

        searchDirection = activeSetDesign - currentDesign;

        nextEqualityLagrangeMultiplier   = activeSetData.EqualityLagrangeMultiplier;
        nextInequalityLagrangeMultiplier = activeSetData.InequalityLagrangeMultiplier(1:numel(inequalityConstraintFunction));

        if isempty(nextEqualityLagrangeMultiplier)
            nextEqualityLagrangeMultiplier = zeros(numel(equalityConstraintFunction), 1);
        end

        % Find the stepsize
        if iteration == 1
            currentWeights.equality = abs(nextEqualityLagrangeMultiplier);
            currentWeights.inequality = abs(nextInequalityLagrangeMultiplier);
        else
            if isempty(equalityConstraintFunction)
                currentWeights.equality = max(abs(nextEqualityLagrangeMultiplier), 0.5 * (previousWeights.equality + abs(nextEqualityLagrangeMultiplier)));
            else
                currentWeights.inequality = max(abs(nextInequalityLagrangeMultiplier), 0.5 * (previousWeights.inequality + abs(nextInequalityLagrangeMultiplier)));
            end
        end


        % Line search with meta function
        [stepSize, stepSizeData, ~] = options.LineSearchFunction( ...
            objectiveFunction, equalityConstraintFunction, inequalityConstraintFunction, ...
            currentDesign, searchDirection, currentWeights);
        
        % Update design variables and enforce bounds
        nextDesign = max(min(currentDesign + stepSize * searchDirection, designSpaceUpperBound), designSpaceLowerBound);

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
            'GradientObjectiveFunction', gradientObjective, ...
            'HessianObjectiveFunction', hessianObjective, ...
            'EqualityLagrangeMultiplier', nextEqualityLagrangeMultiplier, ...
            'InequalityLagrangeMultiplier', nextInequalityLagrangeMultiplier, ...
            'MetaFunctionEqualityWeights', currentWeights.equality,...
            'MetaFunctionInequalityWeights', currentWeights.inequality,...
            'StepSize', stepSize, ...
            'StepSizeData', stepSizeData, ...
            'SearchDirection', searchDirection, ...
            'OptimumCandidatePrevious', currentDesign, ...
            'OptimumCandidateObjectiveValuePrevious', currentObjectiveValue);

        % Update variables for next iteration
        currentDesign = nextDesign;
        currentObjectiveValue = nextObjectiveValue;
        previousWeights = currentWeights;
        currentEqualityLagrangeMultiplier = nextEqualityLagrangeMultiplier;
        currentInequalityLagrangeMultiplier = nextInequalityLagrangeMultiplier;
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
