function [stepSize, stepSizeData] = line_search_golden_ratio(objectiveFunction, startingDesign, searchDirection, maxStepSize, varargin)
%LINE_SEARCH_GOLDEN_RATIO Golden section search with bracketing phase for line search.
%
%   [STEPSIZE, STEPSIZEDATA] = LINE_SEARCH_GOLDEN_RATIO(OBJECTIVEFUNCTION, STARTINGDESIGN,
%   SEARCHDIRECTION, MAXSTEPSIZE, ...) performs a golden section line search
%   with initial bracketing to find the optimal step size along a search direction.
%
%   INPUTS:
%       OBJECTIVEFUNCTION - Function handle to evaluate objective
%       STARTINGDESIGN    - Current design point
%       SEARCHDIRECTION   - Direction for line search
%       MAXSTEPSIZE       - Maximum allowable step size
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'EvaluationObjectiveOptions' - Options for objective evaluation
%       'Alpha'            - Initial step size (default: 0.05)
%       'Beta'             - Expansion factor (default: 1.5)
%       'Rho'              - Golden ratio parameter (default: (1+sqrt(5))/2-1)
%       'Tolerance'        - Convergence tolerance (default: 0.1)
%       'MaxSectioningIterations' - Maximum iterations (default: 20)
%
%   OUTPUTS:
%       STEPSIZE          - Optimal step size found
%       STEPSIZEDATA      - Structure containing search history:
%           .StepSizes    - Array of tested step sizes
%           .ObjectiveValues - Corresponding objective values
%
%   FUNCTIONALITY:
%       1. Bracketing phase to find interval containing minimum
%       2. Golden section search to narrow interval
%       3. Final quadratic approximation for precise step size
%       4. Maintains complete search history for analysis
%
%   NOTES:
%       - Uses golden ratio (0.618) for sectioning
%       - Automatically handles maximum step size constraint
%       - Provides detailed output of search process
%       - Includes polynomial interpolation for final refinement
%
%   REFERENCES:
%       [1] Press, W.H., et al. (2007). "Numerical Recipes", 3rd ed.
%
%   SEE ALSO:
%       line_search_backtracking, evaluate_optimization_objective
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

    p = inputParser;
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Beta', 1.5, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Rho', (1+sqrt(5))/2 - 1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Tolerance', 0.1, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'MaxSectioningIterations', 20, @(x)isnumeric(x)&&isscalar(x));
    parse(p, varargin{:});
    options = p.Results;

    sectioningRatio = min(options.Rho, 1 - options.Rho);
    if(sectioningRatio <= 0)
        error('Sectioning ratio must be positive and less than 1.');
    end

    % **Bracketing Phase**: Find an interval [a, b] where the minimum lies
    stepSizeLeft = 0;
    stepSizeMiddle = 0;
    stepSizeRight = options.Alpha;
    
    % Evaluate objective function at the boundary points
    leftObjective = evaluate_optimization_objective(objectiveFunction, (startingDesign + stepSizeLeft * searchDirection), options.EvaluationObjectiveOptions{:});
    middleObjective = leftObjective;
    rightObjective = evaluate_optimization_objective(objectiveFunction, (startingDesign + stepSizeRight * searchDirection), options.EvaluationObjectiveOptions{:});
    
    % Store step size and objective values in the struct
    stepSizeData.StepSizes = [stepSizeLeft, stepSizeRight];
    stepSizeData.ObjectiveValues = [leftObjective, rightObjective];
    
    hasConverged = false;
    while (~hasConverged)  % Expand the bracket until we find a rise
        stepSizeLeft = stepSizeMiddle;
        leftObjective = middleObjective;

        stepSizeMiddle = stepSizeRight;
        middleObjective = rightObjective;

        stepSizeRight = options.Beta * stepSizeRight;
        if stepSizeRight > maxStepSize  % Prevent growing beyond limit
            stepSizeRight = maxStepSize;
            hasConverged = true;
        end
        rightObjective = evaluate_optimization_objective(objectiveFunction, (startingDesign + stepSizeRight * searchDirection), options.EvaluationObjectiveOptions{:});

        if ((leftObjective >= middleObjective) && (middleObjective <= rightObjective))
            hasConverged = true;
        end
        stepSizeData.StepSizes(end+1) = stepSizeRight;
        stepSizeData.ObjectiveValues(end+1) = rightObjective;
    end

    % **Sectioning Phase**: Golden-section search within [a, b]
    newPoint1 = stepSizeLeft + sectioningRatio * (stepSizeRight - stepSizeLeft);
    newPoint2 = stepSizeRight - sectioningRatio * (stepSizeRight - stepSizeLeft);
    
    objectiveNewPoint1 = evaluate_optimization_objective(objectiveFunction, (startingDesign + newPoint1 * searchDirection), options.EvaluationObjectiveOptions{:});
    objectiveNewPoint2 = evaluate_optimization_objective(objectiveFunction, (startingDesign + newPoint2 * searchDirection), options.EvaluationObjectiveOptions{:});
    
    stepSizeData.StepSizes(end+1) = newPoint1;
    stepSizeData.ObjectiveValues(end+1) = objectiveNewPoint1;
    stepSizeData.StepSizes(end+1) = newPoint2;
    stepSizeData.ObjectiveValues(end+1) = objectiveNewPoint2;

    sectioningIterations = 0;
    while ((stepSizeRight - stepSizeLeft > options.Tolerance) && (sectioningIterations < options.MaxSectioningIterations))  % Continue narrowing interval until tolerance is met
        if objectiveNewPoint1 < objectiveNewPoint2
            stepSizeRight = newPoint2;
            rightObjective = objectiveNewPoint2;

            newPoint2 = newPoint1;
            objectiveNewPoint2 = objectiveNewPoint1;

            newPoint1 = stepSizeLeft + sectioningRatio * (stepSizeRight - stepSizeLeft);
            objectiveNewPoint1 = evaluate_optimization_objective(objectiveFunction, (startingDesign + newPoint1 * searchDirection), options.EvaluationObjectiveOptions{:});
            stepSizeData.StepSizes(end+1) = newPoint1;
            stepSizeData.ObjectiveValues(end+1) = objectiveNewPoint1;
        else
            stepSizeLeft = newPoint1;
            leftObjective = objectiveNewPoint1;

            newPoint1 = newPoint2;
            objectiveNewPoint1 = objectiveNewPoint2;

            newPoint2 = stepSizeRight - sectioningRatio * (stepSizeRight - stepSizeLeft);
            objectiveNewPoint2 = evaluate_optimization_objective(objectiveFunction, (startingDesign + newPoint2 * searchDirection), options.EvaluationObjectiveOptions{:});
            stepSizeData.StepSizes(end+1) = newPoint2;
            stepSizeData.ObjectiveValues(end+1) = objectiveNewPoint2;
        end
        sectioningIterations = sectioningIterations + 1;
    end

    % **Final Step Size**: Polynomial approximation
    if (objectiveNewPoint1 < objectiveNewPoint2)
        stepSizeRight = newPoint2;
        rightObjective = objectiveNewPoint2;

        stepSizeMiddle = newPoint1;
        middleObjective = objectiveNewPoint1;
    else
        stepSizeLeft = newPoint1;
        leftObjective = objectiveNewPoint1;

        stepSizeMiddle = newPoint2;
        middleObjective = objectiveNewPoint2;
    end

    % Polynomial approximation to find the final step size
    polynomialApproximationMatrix = [ones(3,1), [stepSizeLeft; stepSizeMiddle; stepSizeRight], [stepSizeLeft; stepSizeMiddle; stepSizeRight].^2];
    polynomialApproximationCoefficients = polynomialApproximationMatrix \ [leftObjective; middleObjective; rightObjective];
    stepSize = -polynomialApproximationCoefficients(2) / (2 * polynomialApproximationCoefficients(3));

    % check for NaN (e.g., third coefficient is zero)
    if(isnan(stepSize))
        stepSize = 0;
    end
    
    stepSizeData.StepSizes(end+1) = stepSize;
    stepSizeData.ObjectiveValues(end+1) = evaluate_optimization_objective(objectiveFunction, (startingDesign + stepSize * searchDirection), options.EvaluationObjectiveOptions{:});
end
