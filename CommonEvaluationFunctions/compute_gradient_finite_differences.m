function [gradientArray,gradientData] = compute_gradient_finite_differences(functionHandle, designSample, varargin)
%COMPUTE_GRADIENT_FINITE_DIFFERENCES Compute numerical gradient using finite differences
%   This function computes the gradient of an objective function using finite
%   differences for multiple design points simultaneously. It supports various
%   finite difference schemes through specification of evaluation points and
%   step size options.
%
%   [GRADIENTARRAY, GRADIENTDATA] = COMPUTE_GRADIENT_FINITE_DIFFERENCES(OBJECTIVEFUNCTION, DESIGNSAMPLE)
%   returns the gradient(s) and a structure containing data used for gradient computation.
%
%   Output behavior:
%     - If the objective function returns a single output, GRADIENTARRAY is a numeric array of size
%       (nFunctionOutput, nDimension, nSample).
%     - If the objective function returns two outputs, GRADIENTARRAY is a cell array where:
%         GRADIENTARRAY{1} is the gradient for the first output,
%         GRADIENTARRAY{2} is the gradient for the second output.
%       Both have size (nFunctionOutput, nDimension, nSample).
%
%   GRADIENTARRAY = COMPUTE_GRADIENT_FINITE_DIFFERENCES(..., NAME, VALUE, ...) allows specification of
%   additional options through name-value pairs:
%       - 'StepSize'              : double (default: 1e-8). The base step size for finite differences.
%       - 'UseRelativeStepSize'   : logical (default: false). If true, the step size is
%                                     adjusted relative to the magnitude of the design sample.
%       - 'RelativePoints'        : array (default: [-1 0 1]). Stencil points for finite differences.
%       - 'IsObjectiveFunction'   : logical (default: true). If true, the objectiveFunction is treated as a function handle and evaluated as an objective.
%       - 'EvaluationObjectiveOptions' : cell array. Options passed to evaluate_optimization_objective.
%       - 'IsConstraintFunction'  : logical (default: false). If true, the function is treated as a constraint function and evaluated using evaluate_optimization_constraint. When set, the output GRADIENTARRAY is always a cell array with two elements: the first for inequality constraints, the second for equality constraints (even if only one is present).
%       - 'EvaluationConstraintOptions' : cell array. Options passed to evaluate_optimization_constraint when evaluating a constraint function. These options are forwarded to the constraint function and can control its behavior or pass additional parameters.
%
%   Inputs:
%       - OBJECTIVEFUNCTION : function handle
%       - DESIGNSAMPLE : (nSample, nDimension) double
%       - 'StepSize' : double
%       - 'UseRelativeStepSize' : logical
%       - 'RelativePoints' : (1, nPoints) array
%
%   Outputs:
%       - GRADIENTARRAY :
%           * If the function returns one output: (nFunctionOutput, nDimension, nSample) double
%           * If the function returns two outputs: cell array with two elements, each of size (nFunctionOutput, nDimension, nSample)
%       - GRADIENTDATA : struct containing evaluation data and parameters (optional)
%
%   See also compute_objective_function, finite_differences_coefficients, 
%   compute_hessian_finite_differences.
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
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

    % Parse inputs
    p = inputParser;
    addRequired(p, 'functionHandle');
    addRequired(p, 'designSample', @isnumeric);
    addParameter(p, 'StepSize', 1e-8, @isnumeric);
    addParameter(p, 'UseRelativeStepSize', false, @islogical);
    addParameter(p, 'RelativePoints', [-1 0 1], @isnumeric);
    addParameter(p, 'IsObjectiveFunction', false, @islogical);
    addParameter(p, 'EvaluationObjectiveOptions', {});
    addParameter(p, 'IsConstraintFunction', false, @islogical);
    addParameter(p, 'EvaluationConstraintOptions', {});
    parse(p, functionHandle, designSample, varargin{:});
    
    stepSize = p.Results.StepSize;
    relativePoints = p.Results.RelativePoints(:);  % Ensure column vector

    % Get dimensions
    [nSample, nDimension] = size(designSample);

    % adjust step size
    if(isscalar(stepSize))
        stepSize = stepSize*ones(1,nDimension);
    end
    if(p.Results.UseRelativeStepSize)
        stepSize = max(abs(stepSize.*designSample),stepSize);
    else
        stepSize = repmat(stepSize,nSample,1);
    end
    
    % Get finite difference coefficients
    finiteDifferencesCoefficients = finite_differences_coefficients(relativePoints, 1);
    
    % Find non-zero coefficients and their corresponding points
    isNotZeroCoefficient = abs(finiteDifferencesCoefficients) > eps;
    activeCoefficients = finiteDifferencesCoefficients(isNotZeroCoefficient);
    activePoints = relativePoints(isNotZeroCoefficient);
    
    % Find if base point (0) is included and where among active points
    isBasePoint = (activePoints == 0);
    isBasePointIncluded = any(isBasePoint);
    if isBasePointIncluded
        basePointIndex = find(isBasePoint);

        % get perturbation (no 0-stencil) and reorder active coefficients
        perturbationPoints = activePoints([1:basePointIndex-1, basePointIndex+1:end]);
        activeCoefficients = activeCoefficients([basePointIndex,1:basePointIndex-1,basePointIndex+1:end]);
    else
        perturbationPoints = activePoints;
    end
    nPerturbation = length(perturbationPoints);
    
    % Initialize evaluation points matrix (only for non-zero coefficient points)
    totalEvalPoints = nPerturbation*nDimension*nSample;
    if isBasePointIncluded
        totalEvalPoints = totalEvalPoints + nSample;
    end
    evaluationPoints = nan(totalEvalPoints, nDimension);
    
    currentIndex = 1;
    for iDesign = 1:nSample
        % Add base design points if needed
        if isBasePointIncluded
            evaluationPoints(currentIndex,:) = designSample(iDesign,:);
            currentIndex = currentIndex + 1;
        end
        
        % Add perturbation points -> do this procedurally for each variable
        for iDimension = 1:nDimension
            % Number of perturbation points for this variable
            rangeIndex = [currentIndex:(currentIndex + nPerturbation - 1)]';
            
            % Add perturbations only to current variable
            evaluationPoints(rangeIndex,:) = repmat(designSample(iDesign,:),nPerturbation,1);
            evaluationPoints(rangeIndex,iDimension) = evaluationPoints(rangeIndex,iDimension) + stepSize(iDesign,iDimension).*perturbationPoints;
            
            currentIndex = currentIndex + nPerturbation;
        end
    end
    
    % Evaluate objective function at all points
    evaluationData = [];
    secondaryFunctionValues = [];
    isSingleOutput = true;
    if(p.Results.IsObjectiveFunction)
        [primaryFunctionValues,evaluationData] = evaluate_optimization_objective(functionHandle,evaluationPoints,p.Results.EvaluationObjectiveOptions{:});
    else
        % use the constraint evaluation function to deal with any other types of inputs
        % (will also work even if the function is not a constraint)
        [primaryFunctionValues,secondaryFunctionValues,isSingleOutput] = evaluate_optimization_constraint(functionHandle,evaluationPoints,p.Results.EvaluationConstraintOptions{:});
    end

    % if data is not required as output, free up memory 
    isOutputData = (nargout>=2);
    if(~isOutputData)
        evaluationData = [];
    end

    % if the function is a constraint, return as two outputs independently of the number of outputs of the function
    if(p.Results.IsConstraintFunction)
        isSingleOutput = false;
    end

    gradientArray = compute_gradient_finite_differences_given_function_values(primaryFunctionValues, isBasePointIncluded, activeCoefficients, stepSize);
    if(~isSingleOutput)
        gradientArray = {gradientArray,compute_gradient_finite_differences_given_function_values(secondaryFunctionValues, isBasePointIncluded, activeCoefficients, stepSize)};
    end
    
    if(isOutputData)
        gradientData = struct(...
            'RelativePoints',relativePoints,...
            'StepSize',stepSize,...
            'PointEvaluated',evaluationPoints,...
            'FunctionValues',primaryFunctionValues,...
            'SecondaryFunctionValues',secondaryFunctionValues,...
            'ObjectiveEvaluationData',evaluationData,...
            'FiniteDifferencesCoefficients',finiteDifferencesCoefficients);
    end
end