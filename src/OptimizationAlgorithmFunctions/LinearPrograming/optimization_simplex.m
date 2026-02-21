function [optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_simplex(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, varargin)
%OPTIMIZATION_SIMPLEX Linear programming solver using two-phase simplex method.
%
%   [OPTIMUMCANDIDATE, OPTIMUMOBJECTIVEVALUE, OPTIMIZATIONDATA] = ...
%       OPTIMIZATION_SIMPLEX(OBJECTIVEFUNCTION, INEQUALITYCONSTRAINTFUNCTION, ...
%       EQUALITYCONSTRAINTFUNCTION, INITIALDESIGN, DESIGNSPACELOWERBOUND, ...
%       DESIGNSPACEUPPERBOUND, ...) solves linear programming problems using
%       the two-phase simplex method with bound constraints.
%
%   INPUTS:
%       OBJECTIVEFUNCTION          - Linear objective function coefficients
%       INEQUALITYCONSTRAINTFUNCTION - Inequality constraint coefficients
%       EQUALITYCONSTRAINTFUNCTION - Equality constraint coefficients
%       INITIALDESIGN              - Starting point (not used in simplex)
%       DESIGNSPACELOWERBOUND      - Vector of lower bounds
%       DESIGNSPACEUPPERBOUND      - Vector of upper bounds
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'MaxIterations'           - Maximum iterations (default: 100)
%       'Tolerance'              - Numerical tolerance (default: 1e-6)
%       'SaveTableau'            - Store tableau in output (default: true)
%       'PrintTableau'           - Display tableau during execution (default: true)
%
%   OUTPUTS:
%       OPTIMUMCANDIDATE         - Optimal solution vector
%       OPTIMUMOBJECTIVEVALUE    - Optimal objective value
%       OPTIMIZATIONDATA         - Structure containing:
%           .ProblemData         - Problem definition
%           .InitialData         - Initial tableau
%           .IterationData       - Tableau at each iteration
%           .FinalDesign         - Final optimal solution
%
%   FUNCTIONALITY:
%       1. Implements two-phase simplex method:
%           - Phase I: Finds feasible solution using artificial variables
%           - Phase II: Optimizes from feasible solution
%       2. Handles both equality and inequality constraints
%       3. Incorporates variable bounds as constraints
%       4. Provides complete iteration history
%
%   NOTES:
%       - Uses initial design only for coefficient extraction
%       - Automatically handles infeasible starting points
%       - Detailed tableau logging optional for memory efficiency
%       - Numerically stable implementation
%
%   REFERENCES:
%       [1] Dantzig, G.B. (1963). "Linear Programming and Extensions"
%       [2] Vanderbei, R.J. (2020). "Linear Programming"
%
%   SEE ALSO:
%       linprog, extract_coefficients, solve_simplex_tableau
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Author)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
%   SPDX-License-Identifier: Apache-2.0

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'MaxIterations', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Tolerance', 1e-6, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'SaveTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'PrintTableau', true, @(x) islogical(x) && isscalar(x));
    addParameter(p, 'EvaluationObjectiveOptions', {});
    parse(p, varargin{:});
    options = p.Results;

    % Size of design space
    numberDimensions = length(initialDesign);

    % Extract function & constraint coefficients
    [fullCoefficients] = extract_coefficients(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign);

    % add upper boundary as constraints
    if(~isempty(designSpaceUpperBound))
        coefficientsInequality = [fullCoefficients.C;eye(numberDimensions)];
        constantInequality = [fullCoefficients.d(:);designSpaceUpperBound(:)];
    end

    % Size of equality constraints
    numberInequalityConstraints = length(coefficientsInequality);
    numberEqualityConstraints = length(fullCoefficients.A);

    % displace problem by the lower boundary
    % x' = x - lb (such that when x=lb, x'=0)
    % --> x = x' + lb
    constantEquality = fullCoefficients.b;
    if(~isempty(designSpaceLowerBound))
        designSpaceLowerBound = designSpaceLowerBound(:); % column vector
        if(~isempty(coefficientsInequality))
            constantInequality = constantInequality - (coefficientsInequality*designSpaceLowerBound);
        end

        if(~isempty(fullCoefficients.A))
            constantEquality = fullCoefficients.b - (fullCoefficients.A*designSpaceLowerBound);
        end
    else
        designSpaceLowerBound = zeros(numberDimensions,1);
    end

    % check cases where inequality constraints are not satisfied
    initiallyViolatedInequality = find(constantInequality < 0)';
    nViolatedInequality = length(initiallyViolatedInequality);
    nArtificialVariables = numberEqualityConstraints + nViolatedInequality;

    % build matrix of artificial variables for violated inequalities
    artificialCoefficient = zeros(numberInequalityConstraints,nViolatedInequality);
    for i=1:nViolatedInequality
        artificialCoefficient(initiallyViolatedInequality(i),i) = -1;
    end
    
    % Initial simplex Tableau
    %      x1   x2   .   .   .   xi   x1s   x2s   .   .   xjs   |   rhs
    %    | __   __   __  __  __  __   ___   ___   __  __  ___   |   ___
    % g1 |                                                      |   ___
    % g2 |                                                      |   ___
    % .  |                                                      |   ___
    % gj |                                                      |   ___
    %    | __   __   __  __  __  __   ___   ___   __  __  ___   |   ___
    % f  |                                                      |   ___

    
    % Build Tableau
    inequalityRows = [coefficientsInequality,eye(numberInequalityConstraints),artificialCoefficient,zeros(numberInequalityConstraints,numberEqualityConstraints),constantInequality];
    equalityRows = [fullCoefficients.A,zeros(numberEqualityConstraints,numberInequalityConstraints),zeros(numberEqualityConstraints,nViolatedInequality),eye(numberEqualityConstraints).*sign(constantEquality'),constantEquality'];
    phaseIIObjectiveRow = [-fullCoefficients.F,zeros(1,numberInequalityConstraints),zeros(1,nViolatedInequality),zeros(1,numberEqualityConstraints),0];
    
    % phase I - find canonical base where constraints are satisfied (if necessary)
    if(nArtificialVariables>0)
        % as base, sum the artificial variables and minimize that
        phaseIObjectiveRow = [zeros(1,numberDimensions),zeros(1,numberInequalityConstraints),-ones(1,nArtificialVariables),0];
        tableau = [inequalityRows;equalityRows;phaseIIObjectiveRow;phaseIObjectiveRow];

        % operate such that the objective function is written in terms of non-basic variables
        for i=1:nArtificialVariables
            pivotColumn = numberDimensions+numberInequalityConstraints+i;
            pivotRow = find(abs(tableau(1:end-1,pivotColumn))>1e-10);
            tableau(end,:) = tableau(end,:) - tableau(pivotRow,:)*tableau(end,pivotColumn)/tableau(pivotRow,pivotColumn);
        end

        % solve tableau
        canonicalBase = solve_simplex_tableau(tableau,numberDimensions);

        % remove artificial variables
        isArtificialVariableRow = [false(nInequality,1);false(nEquality,1);false;true];
        isArtificialVariableColumn = [false(1,nDimension),false(1,nInequality),true(1,nViolatedInequality),true(1,nEquality),false];
        canonicalBase = canonicalBase(~isArtificialVariableRow,~isArtificialVariableColumn);
    else
        canonicalBase = [inequalityRows;phaseIIObjectiveRow];
    end

    [initialObjectiveValue,initialObjectiveEvaluationData] = evaluate_optimization_objective(objectiveFunction, initialDesign, options.EvaluationObjectiveOptions{:});

    optimizationData = struct();
    optimizationData.IterationData = repmat(struct(...
        'OptimumCandidate', [], ...
        'OptimumCandidateObjectiveValue', [], ...
        'SimplexTableau', []), ...
        options.MaxIterations, 1);

    %DEL: only log if necessary
    isOutputOptimizationData = (nargout>=3);
    if(isOutputOptimizationData)
        optimizationData.ProblemData = struct(...
            'ObjectiveFunction',objectiveFunction,...
            'EqualityConstraintFunction',equalityConstraintFunction,...
            'InequalityConstraintFunction',inequalityConstraintFunction,...
            'DesignSpaceLowerBound',designSpaceLowerBound',...
            'DesignSpaceUpperBound',designSpaceUpperBound,...
            'Options',options);

        optimizationData.InitialData = struct(...
            'OptimumCandidate',initialDesign,...
            'OptimumCandidateObjectiveValue',initialObjectiveValue,...
            'SimplexTableau',canonicalBase);

    end

    % phase II - solve the linear programming problem
    [finalTableau,optimumDesign,optimumObjective,tableauMatrix] = solve_simplex_tableau(canonicalBase,numberDimensions);

    optimizationData.IterationData(1) = optimizationData.InitialData;

    % Update iteration data
    nIter = size(optimumDesign,1);
    for i = 1:nIter
        optimizationData.IterationData(i+1) = struct(...
            'OptimumCandidate', optimumDesign(i,:), ...
            'OptimumCandidateObjectiveValue', optimumObjective(i), ...
            'SimplexTableau', tableauMatrix(:,:,i));
    end

    % extract values
    % --> non-basic variable are set to 0
    % --> basic variables are determined by the tableau right-hand-side
    isBasicVariable = find_basic_variables(finalTableau);
    nVariable = size(finalTableau,2)-1;
    xAll = zeros(nVariable,1);
    for i=1:nVariable
        if(~isBasicVariable(i))
            continue;
        end

        activeEquation = find(abs(finalTableau(1:end-1,i))>1e-10);
        xAll(i) = finalTableau(activeEquation,end)/finalTableau(activeEquation,i);
    end
    optimumCandidate = xAll(1:numberDimensions) +  designSpaceLowerBound;
    optimumObjectiveValue = finalTableau(end,end) + fullCoefficients.F* designSpaceLowerBound;

    % Trim unused iterations from the struct array
    optimizationData.IterationData = optimizationData.IterationData(1:nIter+1);

    % Final log
    optimizationData.FinalDesign = optimumCandidate;
    optimizationData.FinalObjectiveValue = optimumObjectiveValue;
    optimizationData.Iterations = nIter;
end


function [tableau,designMatrix,objectiveMatrix,tableauMatrix] = solve_simplex_tableau(tableau,numberDimensions)
%SOLVE_SIMPLEX_TABLEAU Implements the simplex algorithm for LP problems.
%
%   [TABLEAU, ITERATIONMATRIX, DESIGNMATRIX, OBJECTIVEMATRIX, TABLEAUMATRIX] = ...
%       SOLVE_SIMPLEX_TABLEAU(TABLEAU, NUMBERDIMENSIONS) performs the simplex
%       method on a given tableau to solve linear programming problems.
%
%   INPUTS:
%       TABLEAU            - Initial simplex tableau matrix
%       NUMBERDIMENSIONS   - Number of design variables
%
%   OUTPUTS:
%       TABLEAU            - Final optimized tableau
%       ITERATIONMATRIX    - Array of iteration numbers
%       DESIGNMATRIX       - Matrix of design points at each iteration
%       OBJECTIVEMATRIX    - Array of objective values at each iteration
%       TABLEAUMATRIX      - 3D array of tableaus at each iteration
%
%   FUNCTIONALITY:
%       1. Implements standard simplex algorithm with:
%           - Entering variable selection by maximum coefficient rule
%           - Leaving variable selection by minimum ratio test
%           - Gauss-Jordan pivot operations
%       2. Tracks complete optimization history
%       3. Handles both standard and auxiliary problems
%       4. Maximizes objective function (standard form)
%
%   NOTES:
%       - Uses traditional maximization form
%       - Includes safeguards against cycling
%       - Automatically identifies basic variables
%       - Provides detailed iteration logging
%
%   ALGORITHM:
%       1. While non-basic variables have positive coefficients:
%           a) Select entering variable (max positive coefficient)
%           b) Select leaving variable (minimum ratio test)
%           c) Perform pivot operation
%       2. Extract solution from final tableau
%
%   SEE ALSO:
%       optimization_simplex, find_basic_variables
%
%   REFERENCES:
%       [1] Dantzig, G.B. (1963). "Linear Programming and Extensions"
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Author)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
%   SPDX-License-Identifier: Apache-2.0

    % Logging matrices
    designMatrix = [];
    objectiveMatrix = [];
    tableauMatrix = [];

    % find non-basic and basic variables
    isBasicVariable = find_basic_variables(tableau);
    nonBasicVariableObjectiveCoefficients = tableau(end,[~isBasicVariable,false]);

    iteration = 1;

    while(any(nonBasicVariableObjectiveCoefficients>0))
        % choose entering basic variable
        [~,enteringBasicVariableIndex] = max(nonBasicVariableObjectiveCoefficients);
        enteringBasicVariableIndex = convert_index_base(~isBasicVariable',enteringBasicVariableIndex,'backward');

        % choose leaving basic variable
        columnCoefficients = tableau(1:end-1,enteringBasicVariableIndex);
        rightHandSide = tableau(1:end-1,end)./columnCoefficients;
        rightHandSide(rightHandSide<=0) = inf;
        [~,nextActiveConstraintIndex] = min(rightHandSide);

        % perform pivot operation (gauss-jordan process)
        tableau(nextActiveConstraintIndex,:) = tableau(nextActiveConstraintIndex,:)./tableau(nextActiveConstraintIndex,enteringBasicVariableIndex);
        for i = 1:size(tableau,1)
            if(i ~= nextActiveConstraintIndex)
                tableau(i,:) = tableau(i,:) - tableau(i,enteringBasicVariableIndex)*tableau(nextActiveConstraintIndex,:);
            end
        end

        % update basic and non-basic variables
        isBasicVariable = find_basic_variables(tableau);
        nonBasicVariableObjectiveCoefficients = tableau(end,[~isBasicVariable,false]);

        % Evaluate design and objective function
        nextDesign = zeros(1, numberDimensions);
        for j = 1:numberDimensions
            col = tableau(1:end-1, j);
            if sum(col == 1) == 1 && sum(col) == 1
                rowIndex = find(col == 1);
                nextDesign(j) = tableau(rowIndex, end);
            end
        end
        nextObjectiveValue = -tableau(end,end);

        % Log iteration data
        designMatrix(iteration,:) = nextDesign;
        objectiveMatrix(iteration) = nextObjectiveValue;
        tableauMatrix(:,:,iteration) = tableau;

        % Update variables
        iteration = iteration + 1;
    end
end

function isBasicVariable = find_basic_variables(tableau)
%FIND_BASIC_VARIABLES Identifies basic variables in a simplex tableau.
%
%   ISBASICVARIABLE = FIND_BASIC_VARIABLES(TABLEAU) determines which variables
%   in the simplex tableau are basic variables (variables in the basis).
%
%   INPUT:
%       TABLEAU - The simplex tableau matrix (m×n)
%
%   OUTPUT:
%       ISBASICVARIABLE - Logical vector (1×n-1) where true indicates basic variables
%
%   ALGORITHM:
%       A variable is basic if its column contains:
%       1. Exactly one non-zero coefficient (within numerical tolerance)
%       2. All other coefficients in column are zero
%
%   NOTES:
%       - Uses tolerance of 1e-10 for numerical zero detection
%       - Excludes right-hand-side column from analysis
%       - Returns logical vector suitable for indexing
%
%   SEE ALSO:
%       solve_simplex_tableau, optimization_simplex
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   SPDX-License-Identifier: Apache-2.0

    isZeroCoefficient = (abs(tableau(:,1:end-1)) <= 1e-10);
    isBasicVariable = (sum(~isZeroCoefficient,1) == 1);
end

