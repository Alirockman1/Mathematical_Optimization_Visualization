function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_active_set(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
%OPTIMIZATION_ACTIVE_SET Active set method for constrained optimization.
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE, OPTIMIZATIONDATA] = ...
%       OPTIMIZATION_ACTIVE_SET(OBJECTIVEFUNCTION, INEQUALITYCONSTRAINTFUNCTION, ...
%       EQUALITYCONSTRAINTFUNCTION, INITIALDESIGN, DESIGNSPACELOWERBOUND, ...
%       DESIGNSPACEUPPERBOUND, ...) solves constrained optimization problems
%       using the active set method, handling both equality and inequality
%       constraints with variable bounds.
%
%   INPUTS:
%       OBJECTIVEFUNCTION          - Function handle for objective (or quadratic struct)
%       INEQUALITYCONSTRAINTFUNCTION - Function handle for inequality constraints (or C,d matrices)
%       EQUALITYCONSTRAINTFUNCTION - Function handle for equality constraints (or A,b matrices)
%       INITIALDESIGN              - Starting point vector
%       DESIGNSPACELOWERBOUND      - Vector of lower bounds
%       DESIGNSPACEUPPERBOUND      - Vector of upper bounds
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'MaxIterations'           - Maximum iterations (default: 100)
%       'LagrangeTolerance'      - Tolerance for constraint activity (default: 1e-10)
%       'Tolerance'              - Convergence tolerance (default: 1e-10)
%       'EvaluationObjectiveOptions' - Options for objective evaluation
%       'LineSearchOptions'      - Options for line search
%
%   OUTPUTS:
%       OPTIMUMCANDIDATE         - Optimal solution found
%       OPTIMUMOBJECTIVEVALUE    - Objective value at optimum
%       OPTIMIZATIONDATA         - Structure containing optimization history:
%           .ProblemData         - Problem definition
%           .InitialData         - Initial point information
%           .IterationData       - Detailed iteration records
%           .FinalDesign         - Final optimal solution
%           .LagrangeMultipliers - Final multiplier values
%
%   FUNCTIONALITY:
%       1. Handles both quadratic programming and general nonlinear problems
%       2. Manages active constraints through KKT conditions
%       3. Supports both equality and inequality constraints
%       4. Incorporates design variable bounds
%       5. Provides detailed iteration logging
%
%   NOTES:
%       - Automatically detects quadratic vs. nonlinear problem formulation
%       - Uses Lagrange multipliers to manage constraint activity
%       - Implements constraint dropping when multipliers become negative
%       - Maintains complete optimization history when requested
%
%   REFERENCES:
%       [1] Nocedal, J., Wright, S.J. (2006). "Numerical Optimization"
%
%   SEE ALSO:
%       optimization_conjugate_gradient, optimization_newton
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
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
    addParameter(p, 'LagrangeTolerance', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SaveTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'PrintTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'LineSearchOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    % Number of dimensions
    numberDimensions = length(initialDesign);
    numberInequalityConstraints = numel(inequalityConstraintFunction);

    if ~isempty(equalityConstraintFunction)
        numberEqualityConstraints = length(equalityConstraintFunction(initialDesign));
    else
        numberEqualityConstraints = 0;
    end

    % Detect input type: objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction
    if isa(objectiveFunction, 'function_handle') || isa(inequalityConstraintFunction, 'function_handle') || isa(equalityConstraintFunction, 'function_handle')

        % Extract and display function coefficients
        [fullCoefficients] = lagrange_newton(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign);
    
    else

        fullCoefficients.Q = objectiveFunction.Q;
        fullCoefficients.e = objectiveFunction.e;
        fullCoefficients.A = equalityConstraintFunction.A;
        fullCoefficients.b = equalityConstraintFunction.b;
        fullCoefficients.C = inequalityConstraintFunction.C;
        fullCoefficients.d = inequalityConstraintFunction.d;

    end

    C_bound = [];
    d_bound = [];
    if ~isempty(designSpaceLowerBound)
            C_bound = [C_bound; -eye(numberDimensions)];   % lb:  x >= lb  -->  -x <= -lb
            d_bound = [d_bound; -designSpaceLowerBound(:)];
    end
    if ~isempty(designSpaceUpperBound)
        C_bound = [C_bound; eye(numberDimensions)];    % ub:  x <= ub
        d_bound = [d_bound; designSpaceUpperBound(:)];
    end


    % Combine the inequality constraints from the function with the bounds.
    if isempty(fullCoefficients.C)
        coefficientsInequality = C_bound;
        constantInequality = d_bound;
    else
        coefficientsInequality = [fullCoefficients.C; C_bound];
        constantInequality     = [fullCoefficients.d; d_bound];
    end


    inequalityActive = [];

    % For inequality constraints, add those that are active.
    if ~isempty(coefficientsInequality)
        inequalityActive = find(-coefficientsInequality*initialDesign' + constantInequality <= options.LagrangeTolerance);
    end

    % Evaluate the objective function
    [initialObjectiveValue,initialObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, initialDesign, options.EvaluationObjectiveOptions{:});

    optimizationData = struct();
    optimizationData.IterationData = repmat(struct(...
        'Iteration', [], ...
        'OptimumCandidate', [], ...
        'OptimumCandidateObjectiveValue', [], ...
        'OptimumCandidateObjectiveEvaluationData', [], ...
        'InequalityLagrangeMultiplier', [], ...
        'EqualityLagrangeMultiplier', [], ...
        'StepSize', [], ...
        'StepSizeData', [], ...
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

    iteration = 1;
    hasConverged = false;
    currentDesign = initialDesign;

    while ~hasConverged && iteration <= options.MaxIterations
        
        [currentObjectiveValue,currentObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, currentDesign, options.EvaluationObjectiveOptions{:});

        % Determining active inequality constraints
        activeInequalityCoefficients = [];
        if ~isempty(inequalityActive)
            activeInequalityCoefficients = coefficientsInequality(inequalityActive,:);
        end
        numberActiveInequalityConstraints = size(activeInequalityCoefficients,1);

        % Step 0: Assemble the lagrange matrix
        grad = fullCoefficients.Q*currentDesign' + fullCoefficients.e;

        if ~isempty(equalityConstraintFunction)
            LagrangianMatrix = [fullCoefficients.Q, fullCoefficients.A', activeInequalityCoefficients'; 
                                fullCoefficients.A, zeros(size(fullCoefficients.A,1)), zeros(size(fullCoefficients.A,1), size(activeInequalityCoefficients,1)); 
                                activeInequalityCoefficients, zeros(size(activeInequalityCoefficients,1), size(fullCoefficients.A,1)), zeros(size(activeInequalityCoefficients,1))];
            rhs = [-grad; zeros(size(fullCoefficients.A,1)); zeros(size(activeInequalityCoefficients,1))];
        else
            LagrangianMatrix = [fullCoefficients.Q, activeInequalityCoefficients'; 
                                activeInequalityCoefficients, zeros(size(activeInequalityCoefficients,1))];
            rhs = [-grad; zeros(size(activeInequalityCoefficients,1), 1)];
        end

        % Step 1: Compute lagrange multipliers
        lagrangeVector = LagrangianMatrix\rhs;

        % Step 2: Compute the search directions based on the psudodesign
        searchDirectionIndex = 1:numberDimensions;
        equalityLagrangeIndex = (numberDimensions+1):(numberDimensions+numberEqualityConstraints);
        inequalityLagrangeIndex = (numberDimensions+numberEqualityConstraints+1):(numberDimensions+numberEqualityConstraints+numberActiveInequalityConstraints);

        searchDirection = lagrangeVector(searchDirectionIndex,1)';
        equalityLagrangeMultiplier = lagrangeVector(equalityLagrangeIndex,1);
        inequalityLagrangeMultiplier = lagrangeVector(inequalityLagrangeIndex,1);

        % Step 3: Check KKT condition
        if norm(searchDirection) < options.Tolerance
            
            % If any multiplier for an inequality constraint is negative,
            % remove the one with the most negative multiplier.
            if ~isempty(inequalityLagrangeMultiplier) && any(inequalityLagrangeMultiplier < -options.Tolerance)
                [~, idx] = min(inequalityLagrangeMultiplier);
                inequalityActive(idx) = [];                                % Remove the blocking constraint.
                continue;
            else
                % Optimality conditions are satisfied.
                nextDesign = currentDesign;
                nextObjectiveValue = currentObjectiveValue;
                nextObjectiveEvaluationData = currentObjectiveEvaluationData;
                hasConverged = true;
                break;
            end

        else
            % Compute the maximum step length alpha such that no inactive
            % inequality constraint is violated.
            maxStepSize = 1;
            blockingConstraint = [];
            stepSizeArray = [];
            if ~isempty(coefficientsInequality)
                % Consider constraints not in the working set.
                nonActive = setdiff((1:size(coefficientsInequality,1))', inequalityActive);
                for i = nonActive'
                    % Only consider constraints that would be violated when moving along p.
                    if coefficientsInequality(i,:) * (currentDesign +...
                            maxStepSize*searchDirection)' -...
                            constantInequality(i) > options.Tolerance
                        iStepSize = (constantInequality(i) -...
                            coefficientsInequality(i,:)*currentDesign') /...
                            (coefficientsInequality(i,:)*searchDirection');
                        if iStepSize < maxStepSize
                            maxStepSize = iStepSize;
                            blockingConstraint = i;
                        end
                    end
                    stepSizeArray = [stepSizeArray,maxStepSize];
                end
            end
            
            % Update design
            nextDesign = currentDesign + maxStepSize*searchDirection;

            % If a blocking constraint was encountered (step < full step),
            % add it to the working set.
            if maxStepSize < 1 - options.Tolerance && ~isempty(blockingConstraint)
                inequalityActive = [inequalityActive; blockingConstraint];
            end
        end

        [nextObjectiveValue,nextObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, nextDesign, options.EvaluationObjectiveOptions{:});


        % Log iteration data
        optimizationData.IterationData(iteration) = struct(...
            'Iteration', iteration, ...
            'OptimumCandidate', nextDesign, ...
            'OptimumCandidateObjectiveValue', nextObjectiveValue, ...
            'OptimumCandidateObjectiveEvaluationData', nextObjectiveEvaluationData, ...
            'InequalityLagrangeMultiplier', inequalityLagrangeMultiplier, ...
            'EqualityLagrangeMultiplier', equalityLagrangeMultiplier, ...
            'StepSize', maxStepSize, ...
            'StepSizeData', stepSizeArray, ...
            'SearchDirection', searchDirection, ...
            'OptimumCandidatePrevious', currentDesign, ...
            'OptimumCandidateObjectiveValuePrevious', currentObjectiveValue);

        currentDesign = nextDesign;
        
        iteration = iteration +1;

    end

    % Trim unused iterations from the struct array
    optimizationData.IterationData = optimizationData.IterationData(1:iteration - 1);

    % Outputs
    optimumCandidate = nextDesign;
    optimumObjectiveValue = nextObjectiveValue;
    optimumInequalityLagrangeMultiplier = [inequalityLagrangeMultiplier; ...
    zeros(size(coefficientsInequality, 1) - length(inequalityLagrangeMultiplier), 1)];

    % Final log
    optimizationData.FinalDesign = optimumCandidate;
    optimizationData.EqualityLagrangeMultiplier = equalityLagrangeMultiplier;
    optimizationData.InequalityLagrangeMultiplier = optimumInequalityLagrangeMultiplier;
    optimizationData.FinalObjectiveValue = optimumObjectiveValue;
    optimizationData.FinalObjectiveValueData = nextObjectiveEvaluationData;
    optimizationData.Iterations = iteration - 1;
end

function [fullCoefficients] = lagrange_newton(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign)
%LAGRANGE_NEWTON Extracts quadratic programming coefficients from nonlinear functions.
%
%   FULLCOEFFICIENTS = LAGRANGE_NEWTON(OBJECTIVEFUNCTION, INEQUALITYCONSTRAINTFUNCTION,
%   EQUALITYCONSTRAINTFUNCTION, INITIALDESIGN) converts nonlinear optimization
%   problem components into quadratic programming (QP) formulation coefficients.
%
%   INPUTS:
%       OBJECTIVEFUNCTION          - Function handle for objective function
%       INEQUALITYCONSTRAINTFUNCTION - Function handle/cell array for inequality constraints
%       EQUALITYCONSTRAINTFUNCTION - Function handle/cell array for equality constraints
%       INITIALDESIGN              - Initial design point (determines variable count)
%
%   OUTPUTS:
%       FULLCOEFFICIENTS - Structure containing QP coefficients:
%           .Q          - Quadratic term matrix of objective
%           .e          - Linear term vector of objective  
%           .A          - Equality constraint coefficient matrix
%           .b          - Equality constraint constant vector
%           .C          - Inequality constraint coefficient matrix
%           .d          - Inequality constraint constant vector
%
%   FUNCTIONALITY:
%       1. Uses symbolic differentiation to extract quadratic coefficients
%       2. Handles both single and multiple constraints (cell arrays)
%       3. Converts nonlinear constraints to linear QP form
%       4. Maintains proper sign conventions for QP formulation
%
%   NOTES:
%       - Requires Symbolic Math Toolbox
%       - Automatically handles empty equality constraints
%       - Preserves constraint relationships through sign adjustment
%       - Returns matrices in standard QP format (min 0.5x'Qx + e'x)
%
%   SEE ALSO:
%       optimization_active_set, hessian, jacobian
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   SPDX-License-Identifier: Apache-2.0

    syms x [1 length(initialDesign)] real;

    nDesign = length(initialDesign);

    objExpr = objectiveFunction(x);

    if iscell(inequalityConstraintFunction)
        inequalityConstraintsExpr = sym(zeros(length(inequalityConstraintFunction), 1));

        for i = 1:length(inequalityConstraintFunction)
            inequalityConstraintsExpr(i) = inequalityConstraintFunction{i}(x);
        end

    else
        inequalityConstraintsExpr = inequalityConstraintFunction(x);
    end

    hasEqualityConstraints = ~isempty(equalityConstraintFunction);
    if hasEqualityConstraints
            if iscell(equalityConstraintFunction)
                equalityConstraintsExpr = sym(zeros(length(equalityConstraintFunction), 1));
        
                for i = 1:length(equalityConstraintFunction)
                    equalityConstraintsExpr(i) = equalityConstraintFunction{i}(x);
                end
        
            else
                equalityConstraintsExpr = equalityConstraintFunction(x);
            end
    end
    
    % Create Q matrix of objective function -> Hessian
    coefficientsObjective = double(hessian(objExpr, x));

    % Create e vector of objective function -> Hessian
    constantObjective = double(subs(jacobian(objExpr, x), x, zeros(1, nDesign))');

    % Create inequality matrix and inequality constant vector -> C & d
    % Initialize coefficient matrix with zeros
    coefficientsInequality = zeros(length(inequalityConstraintsExpr), length(x));
    constantInequality = zeros(length(inequalityConstraintsExpr), 1);
    
    % fprintf('\nConstraint Coefficients:\n');
    for i = 1:length(inequalityConstraintsExpr)
        % Extract coefficients of x1 and x2
        [c_terms, c_vars] = coeffs(inequalityConstraintsExpr(i), x);
        
        % Assign coefficients to appropriate positions
        for j = 1:length(c_vars)
            varIndex = find(x == c_vars(j), 1);                            % Find the variable index
            if ~isempty(varIndex)
                coefficientsInequality(i, varIndex) = c_terms(j);          % Assign coefficient
            end
        end
        % Extract constant term (independent term)
        constantTerm = subs(inequalityConstraintsExpr(i), x, [0, 0]); 
        constantInequality(i) = double(constantTerm); % Store in b vector
    end

    % Create equality matrix and equality constant vector -> A & b
    if hasEqualityConstraints
        numEqConstraints = length(equalityConstraintsExpr);
        coeffsEqCon = zeros(numEqConstraints, length(x));
        constantEquality = zeros(numEqConstraints, 1);

        % fprintf('\nExtracting Equality Constraints:\n');
        for i = 1:numEqConstraints
            [c_terms, c_vars] = coeffs(equalityConstraintsExpr(i), x);

            for j = 1:length(c_vars)
                varIndex = find(x == c_vars(j), 1);                        % Find variable index
                if ~isempty(varIndex)
                    coeffsEqCon(i, varIndex) = c_terms(j);
                end
            end

            % Extract constant term
            constantTerm = subs(equalityConstraintsExpr(i), x, zeros(1, length(initialDesign)));
            constantEquality(i) = double(constantTerm);
        end

        coefficientsEquality = coeffsEqCon;
    else
        coefficientsEquality = [];
        constantEquality     = [];
    end
    
    % Store Outputs in Structs
    fullCoefficients = struct('Q', coefficientsObjective, 'e', constantObjective, 'A', -1*coefficientsEquality, 'b', constantEquality, 'C', coefficientsInequality, 'd', -1*constantInequality);

end

