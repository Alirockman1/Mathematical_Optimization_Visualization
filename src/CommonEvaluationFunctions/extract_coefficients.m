function [fullCoefficients] = extract_coefficients(objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, varargin)
%EXTRACT_COEFFICIENTS Extracts coefficients from optimization problem components.
%
%   FULLCOEFFICIENTS = EXTRACT_COEFFICIENTS(OBJECTIVEFUNCTION, ...
%   INEQUALITYCONSTRAINTFUNCTION, EQUALITYCONSTRAINTFUNCTION, INITIALDESIGN, ...)
%   converts symbolic optimization problem components into numerical coefficient
%   matrices for linear or quadratic programming formulations.
%
%   INPUTS:
%       OBJECTIVEFUNCTION          - Function handle for objective
%       INEQUALITYCONSTRAINTFUNCTION - Function handle/cell array for inequalities
%       EQUALITYCONSTRAINTFUNCTION - Function handle/cell array for equalities
%       INITIALDESIGN              - Design point (determines variable count)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'Programming'            - Problem type ('linear' or 'quadratic')
%
%   OUTPUTS:
%       FULLCOEFFICIENTS         - Structure containing:
%           .Q                  - Quadratic term matrix (QP only)
%           .F                  - Linear coefficient vector
%           .e                  - Constant term
%           .A                  - Equality constraint matrix
%           .b                  - Equality constraint vector
%           .C                  - Inequality constraint matrix
%           .d                  - Inequality constraint vector
%
%   FUNCTIONALITY:
%       1. Handles both linear and quadratic programming problems
%       2. Supports function handles and cell arrays of constraints
%       3. Uses symbolic differentiation for coefficient extraction
%       4. Maintains standard forms (Ax=b, Cx≤d)
%       5. Automatically handles empty constraints
%
%   NOTES:
%       - Requires Symbolic Math Toolbox
%       - For QP: Extracts Hessian (Q) and gradient (F) terms
%       - For LP: Returns empty Q matrix
%       - Automatically converts constraints to standard form
%
%   SEE ALSO:
%       sym, hessian, jacobian, optimization_simplex
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

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'Programming', 'linear', @(x) any(validatestring(x, {'linear', 'quadratic'})));
    parse(p, varargin{:});
    options = p.Results;

    % Define symbolic variables
    n = length(initialDesign);
    x = sym('x', [1 n], 'real');

    %% Objective Function
    objExpr = expand(objectiveFunction(x));

    if strcmpi(options.Programming, 'quadratic')
        % Q: Hessian
        quadraticCoefficientsObjective = double(hessian(objExpr, x));
        % e: gradient
        coefficientsObjective = double(subs(jacobian(objExpr, x), x, zeros(1, n)));
        constantObjective = double(subs(objExpr, x, zeros(1, n)));
    else
        % For linear program: No Hessian
        quadraticCoefficientsObjective = [];
        coefficientsObjective = zeros(1, n);
        for j = 1:n
            coefficientsObjective(j) = double(subs(diff(objExpr, x(j)), x, zeros(1, n)));
        end
        constantObjective = double(subs(objExpr, x, zeros(1, n)));
    end

    %% Inequality Constraints (C * x <= d)
    if isempty(inequalityConstraintFunction)
        coefficientsInequality = [];
        constantInequality = [];
    else
        if iscell(inequalityConstraintFunction)
            ineqExprs = [];
            for i = 1:length(inequalityConstraintFunction)
                ineqExprs = [ineqExprs; inequalityConstraintFunction{i}(x)];
            end
        else
            ineqExprs = inequalityConstraintFunction(x);
        end
        numIneq = length(ineqExprs);
        coefficientsInequality = zeros(numIneq, n);
        constantInequality = zeros(numIneq, 1);
        for i = 1:numIneq
            expr = expand(ineqExprs(i));
            for j = 1:n
                coefficientsInequality(i, j) = double(subs(diff(expr, x(j)), x, zeros(1, n)));
            end
            constantInequality(i) = double(subs(expr, x, zeros(1, n)));
        end
    end

    %% Equality Constraints (A * x = b)
    if isempty(equalityConstraintFunction)
        coefficientsEquality = [];
        constantEquality = [];
    else
        if iscell(equalityConstraintFunction)
            equation = [];
            for i = 1:length(equalityConstraintFunction)
                equation = [equation; equalityConstraintFunction{i}(x)];
            end
        else
            equation = equalityConstraintFunction(x);
        end
        numEq = length(equation);
        coefficientsEquality = zeros(numEq, n);
        constantEquality = zeros(numEq, 1);
        for i = 1:numEq
            expr = expand(equation(i));
            for j = 1:n
                coefficientsEquality(i, j) = double(subs(diff(expr, x(j)), x, zeros(1, n)));
            end
            if strcmpi(programType, 'linear')
                constantEquality(i) = 0;  % For linear LP case: b = 0 (homogeneous form)
            else
                constantEquality(i) = double(subs(expr, x, zeros(1, n)));
            end
        end
    end

    % Store Outputs in Structs
    fullCoefficients = struct('Q', quadraticCoefficientsObjective,'F', coefficientsObjective, 'e', constantObjective, 'A', coefficientsEquality, 'b', constantEquality*-1, 'C', coefficientsInequality, 'd', constantInequality*-1);

end
