function [hessianArray,hessianData,gradientArray] = compute_hessian_finite_differences(functionHandle, designSample, varargin)
%COMPUTE_HESSIAN Compute numerical Hessian using finite differences
%   This function computes the Hessian matrix of an objective function using 
%   finite differences for multiple design points simultaneously. It supports 
%   various finite difference schemes through specification of evaluation points
%   and step size options.
%
%   HESSIANARRAY = COMPUTE_HESSIAN(OBJECTIVEFUNCTION, DESIGNSAMPLE) computes 
%   the Hessian using 3-point central differences with default step size of 1e-4 and 
%   relative points [-1 0 1]. The Hessian array has dimensions 
%   (nDimension, nDimension, nFunctionOutput, nSample). Optionally, you can enable 
%   relative step sizing by setting 'UseRelativeStepSize' to true.
%
%   [HESSIANARRAY, HESSIANDATA, GRADIENTARRAY] = COMPUTE_HESSIAN(OBJECTIVEFUNCTION, DESIGNSAMPLE) 
%   returns the Hessian(s), a structure containing data used for Hessian computation, and optionally the gradients.
%
%   Output behavior:
%     - If the objective function returns a single output, HESSIANARRAY is a numeric array of size
%       (nDimension, nDimension, nFunctionOutput, nSample).
%     - If the objective function returns two outputs, HESSIANARRAY is a cell array where:
%         HESSIANARRAY{1} is the Hessian for the first output,
%         HESSIANARRAY{2} is the Hessian for the second output.
%       Both have size (nDimension, nDimension, nFunctionOutput, nSample).
%
%   HESSIANARRAY = COMPUTE_HESSIAN(...,NAME,VALUE,...) allows specification of
%   additional options through name-value pairs:
%       - 'StepSize'              : double (default: 1e-4). The base step size for finite differences.
%       - 'UseRelativeStepSize'   : logical (default: false). If true, the step size is 
%                                     adjusted relative to the magnitude of the design sample.
%       - 'RelativePoints'        : array (default: [-1 0 1]). Stencil points for finite differences.
%       - 'IsObjectiveFunction'   : logical (default: true). If true, the objectiveFunction is treated as a function handle and evaluated as an objective.
%       - 'EvaluationObjectiveOptions' : cell array. Options passed to evaluate_optimization_objective.
%       - 'IsConstraintFunction'  : logical (default: false). If true, the function is treated as a constraint function and evaluated using evaluate_optimization_constraint. When set, the output HESSIANARRAY is always a cell array with two elements: the first for inequality constraints, the second for equality constraints (even if only one is present).
%       - 'EvaluationConstraintOptions' : cell array. Options passed to evaluate_optimization_constraint when evaluating a constraint function. These options are forwarded to the constraint function and can control its behavior or pass additional parameters.
%
%   Inputs:
%       - OBJECTIVEFUNCTION : function handle
%       - DESIGNSAMPLE : (nSample, nDimension) double
%       - 'StepSize' : double
%       - 'UseRelativeStepSize' : logical
%       - 'RelativePoints' : (1,nPoints) array
%
%   Outputs:
%       - HESSIANARRAY :
%           * If the function returns one output: (nDimension, nDimension, nFunctionOutput, nSample) double
%           * If the function returns two outputs: cell array with two elements, each of size (nDimension, nDimension, nFunctionOutput, nSample)
%       - HESSIANDATA : struct containing evaluation data and parameters (optional)
%       - GRADIENTARRAY :
%           * If the function returns one output: (nFunctionOutput, nDimension, nSample) double (optional)
%           * If the function returns two outputs: cell array with two elements, each of size (nFunctionOutput, nDimension, nSample) (optional)
%
%   See also compute_objective_function, finite_differences_coefficients, 
%   compute_gradient_finite_differences, compute_hessian_given_function_values.
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
    addParameter(p, 'StepSize', 1e-4, @isnumeric);
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
    
    % Get finite difference coefficients for first and second derivatives
    firstDerivativeCoefficients = finite_differences_coefficients(relativePoints, 1);
    secondDerivativeCoefficients = finite_differences_coefficients(relativePoints, 2);
    
    % Find non-zero coefficients and their corresponding points
    isZeroCoefficient = ((abs(firstDerivativeCoefficients) <= eps) & (abs(secondDerivativeCoefficients) <= eps));
    activeFirstCoefficients = firstDerivativeCoefficients(~isZeroCoefficient);
    activeFirstPoints = relativePoints(~isZeroCoefficient);
    activeSecondCoefficients = secondDerivativeCoefficients(~isZeroCoefficient);
    activeSecondPoints = relativePoints(~isZeroCoefficient);
    
    % Find if base point (0) is included and where among active points
    isBasePointFirst = (activeFirstPoints == 0);
    if(any(isBasePointFirst))
        basePointFirstIndex = find(isBasePointFirst);

        perturbationPointsFirst = activeFirstPoints([1:basePointFirstIndex-1, basePointFirstIndex+1:end]);
        activeFirstCoefficients = activeFirstCoefficients([basePointFirstIndex,1:basePointFirstIndex-1,basePointFirstIndex+1:end]);
    else
        perturbationPointsFirst = activeFirstPoints;
    end
    nPerturbationFirst = length(perturbationPointsFirst);
    
    isBasePointSecond = (activeSecondPoints == 0);
    if(any(isBasePointSecond))
        basePointSecondIndex = find(isBasePointSecond);

        perturbationPointsSecond = activeSecondPoints([1:basePointSecondIndex-1, basePointSecondIndex+1:end]);
        activeSecondCoefficients = activeSecondCoefficients([basePointSecondIndex,1:basePointSecondIndex-1,basePointSecondIndex+1:end]);
    else
        perturbationPointsSecond = activeSecondPoints;
    end
    nPerturbationSecond = length(perturbationPointsSecond);
    
    % Initialize evaluation points matrix
    % We need points for:
    % 1. Base points (if included)
    % 2. Points for diagonal terms (second derivatives)
    % 3. Points for off-diagonal terms (mixed derivatives using tensor product)
    
    % Count total evaluation points needed
    totalEvaluationPoints = 0;
    
    % Base points (if needed)
    isBasePointIncluded = any(isBasePointSecond) || any(isBasePointFirst);
    if(isBasePointIncluded)
        totalEvaluationPoints = totalEvaluationPoints + nSample;
    end
    
    % Diagonal terms (second derivatives)
    totalEvaluationPoints = totalEvaluationPoints + nSample * nDimension * nPerturbationSecond;
    
    % Off-diagonal terms (mixed derivatives)
    % For each pair of dimensions, we need a tensor product of perturbation points
    nPairs = nDimension*(nDimension-1)/2;  % Number of unique dimension pairs
    totalEvaluationPoints = totalEvaluationPoints + nSample*nPairs*(nPerturbationFirst^2);
    
    % Initialize evaluation points matrix
    evaluationPoints = zeros(totalEvaluationPoints, nDimension);
    
    currentIndex = 1;
    for iDesign = 1:nSample
        % Add base design points if needed
        if(isBasePointIncluded)
            evaluationPoints(currentIndex,:) = designSample(iDesign,:);
            currentIndex = currentIndex + 1;
        end

        % Add perturbation points -> do this procedurally for each variable
        for iDimension = 1:nDimension
            % Number of perturbation points for this variable
            rangeIndex = [currentIndex:(currentIndex + nPerturbationSecond - 1)]';
            
            % Add perturbations only to current variable
            evaluationPoints(rangeIndex,:) = repmat(designSample(iDesign,:),nPerturbationSecond,1);
            evaluationPoints(rangeIndex,iDimension) = evaluationPoints(rangeIndex,iDimension) + stepSize(iDesign,iDimension).*perturbationPointsSecond;
            
            currentIndex = currentIndex + nPerturbationSecond;
        end
        
        % Add points for off-diagonal terms (mixed derivatives)
        % For each unique pair of dimensions
        for iVariable = 1:nDimension
            for jVariable = (iVariable+1):nDimension
                % Create tensor product of perturbation points
                for iPerturbation = 1:nPerturbationFirst
                    for jPerturbation = 1:nPerturbationFirst
                        % Add perturbations to both variables
                        evaluationPoints(currentIndex,iVariable) = designSample(iDesign,iVariable) + stepSize(iDesign,iVariable) * perturbationPointsFirst(iPerturbation);
                        evaluationPoints(currentIndex,jVariable) = designSample(iDesign,jVariable) + stepSize(iDesign,jVariable) * perturbationPointsFirst(jPerturbation);
                        
                        currentIndex = currentIndex + 1;
                    end
                end
            end
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

    hessianArray = compute_hessian_finite_differences_given_function_values(primaryFunctionValues, isBasePointIncluded, activeFirstCoefficients, activeSecondCoefficients, stepSize);
    if(~isSingleOutput)
        hessianArray = {hessianArray,compute_hessian_finite_differences_given_function_values(secondaryFunctionValues, isBasePointIncluded, activeFirstCoefficients, activeSecondCoefficients, stepSize)};
    end

    if(isOutputData)
        hessianData = struct(...
            'RelativePoints', relativePoints, ...
            'StepSize', stepSize, ...
            'PointsEvaluated', evaluationPoints, ...
            'FunctionValues', primaryFunctionValues, ...
            'SecondaryFunctionValues', secondaryFunctionValues, ...
            'ObjectiveEvaluationData', evaluationData, ...
            'FirstDerivativeCoefficients', activeFirstCoefficients, ...
            'SecondDerivativeCoefficients', activeSecondCoefficients);
    end

    % use already computed data to estimate gradient as well
    if(nargout>=3)
        firstDerivativeRelevant = exclude_mixed_derivative_terms(isBasePointIncluded, nSample, nDimension, nPerturbationFirst, nPerturbationSecond);
        
        firstDerivativeValuesPrimary = get_array_subset_if_nonempty(primaryFunctionValues,firstDerivativeRelevant);
        gradientArray = compute_gradient_finite_differences_given_function_values(firstDerivativeValuesPrimary, isBasePointIncluded, activeFirstCoefficients, stepSize);

        if(~isSingleOutput)
            firstDerivativeValuesSecondary = get_array_subset_if_nonempty(secondaryFunctionValues,firstDerivativeRelevant);
            gradientArray = {gradientArray,compute_gradient_finite_differences_given_function_values(firstDerivativeValuesSecondary, isBasePointIncluded, activeFirstCoefficients, stepSize)};
        end
    end
end

function isFirstDerivativeRelevant = exclude_mixed_derivative_terms(isBasePointIncluded, nSample, nDimension, nPerturbationFirst, nPerturbationSecond)
    isFirstDerivativeRelevant = false(nSample,1);
    
    blockSizeIndividual = nPerturbationSecond*nDimension;
    blockSizeMixed = nPerturbationFirst^2*nDimension*(nDimension-1)/2;
    currentIndex = 1;
    for iSample = 1:nSample
        if(isBasePointIncluded)
            isFirstDerivativeRelevant(currentIndex) = true;
            currentIndex = currentIndex + 1;
        end
        
        rangeIndex = [currentIndex:(currentIndex + blockSizeIndividual - 1)]';
        isFirstDerivativeRelevant(rangeIndex) = true;
        currentIndex = currentIndex + blockSizeIndividual;

        rangeIndex = [currentIndex:(currentIndex + blockSizeMixed - 1)]';
        isFirstDerivativeRelevant(rangeIndex) = false;
        currentIndex = currentIndex + blockSizeMixed;
    end
end

function hessianArray = compute_hessian_finite_differences_given_function_values(functionValues, isBasePointIncluded, activeCoefficientsFirst, activeCoefficientsSecond, stepSize)
%COMPUTE_HESSIAN_GIVEN_FUNCTION_VALUES Compute Hessian from pre-evaluated function values
%   This function computes the Hessian matrix using pre-evaluated function values at
%   perturbation points. It is primarily used as a helper function by other
%   optimization routines that need to compute Hessians efficiently from
%   existing function evaluations.
%
%   HESSIANARRAY = COMPUTE_HESSIAN_GIVEN_FUNCTION_VALUES(FUNCTIONVALUES,
%   ISBASEPOINTINCLUDED, ACTIVECOEFFICIENTS_FIRST, ACTIVECOEFFICIENTS_SECOND,
%   STEPSIZE) computes the Hessian using pre-evaluated function values and
%   finite difference coefficients for both first and second derivatives.
%
%   Inputs:
%       - FUNCTIONVALUES : (nPoints, nFunctionOutput) double
%           Pre-evaluated function values at perturbation points
%       - ISBASEPOINTINCLUDED : logical
%           Whether the base point (unperturbed) evaluation is included
%       - ACTIVECOEFFICIENTS_FIRST : (nCoefficients, 1) double
%           Finite difference coefficients for first derivatives (mixed terms)
%       - ACTIVECOEFFICIENTS_SECOND : (nCoefficients, 1) double
%           Finite difference coefficients for second derivatives (diagonal terms)
%       - STEPSIZE : (nSample, nDimension) double
%           Step sizes used for perturbations in each dimension
%
%   Output:
%       - HESSIANARRAY : (nDimension, nDimension, nFunctionOutput, nSample) double
%           Computed Hessian array containing second derivatives
%
%   See also compute_hessian_finite_differences, compute_gradient_given_function_values,
%   finite_differences_coefficients.
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

    if(isempty(functionValues))
        hessianArray = [];
        return;
    end

    % Extract and reshape values
    nFunctionOutput = size(functionValues,2);
    nSample = size(stepSize,1);
    nDimension = size(stepSize,2);
    
    nPerturbationFirst = length(activeCoefficientsFirst);
    if(isBasePointIncluded)
        nPerturbationFirst = nPerturbationFirst - 1;
    end

    nPerturbationSecond = length(activeCoefficientsSecond);
    if(isBasePointIncluded)
        nPerturbationSecond = nPerturbationSecond - 1;
    end

    % Initialize Hessian array
    hessianArray = zeros(nDimension, nDimension, nFunctionOutput, nSample);
    blockSize = nPerturbationSecond*nDimension + nPerturbationFirst^2*nDimension*(nDimension-1)/2;
    if(isBasePointIncluded)
        blockSize = blockSize + 1;
    end

    currentIndex = 1;
    for iDesign = 1:nSample
        rangeIndex = [currentIndex:(currentIndex + blockSize - 1)]';
        currentFunctionValues = functionValues(rangeIndex,:);

        for iOutput = 1:nFunctionOutput
            hessianArray(:,:,iOutput,iDesign) = compute_hessian_current_point(currentFunctionValues(:,iOutput), isBasePointIncluded, activeCoefficientsFirst, activeCoefficientsSecond, stepSize(iDesign,:));
        end

        currentIndex = currentIndex + blockSize;
    end
end

function hessianArray = compute_hessian_current_point(functionValues, isBasePointIncluded, activeCoefficientsFirst, activeCoefficientsSecond, stepSize)
    nDimension = size(stepSize,2);
    hessianArray = nan(nDimension, nDimension);

    if(isBasePointIncluded)
        nPerturbationSecond = length(activeCoefficientsSecond) - 1;
        currentIndex = 2;
    else
        nPerturbationSecond = length(activeCoefficientsSecond);
        currentIndex = 1;
    end

    % Compute diagonal terms (second derivatives)
    currentFunctionValues = nan(length(activeCoefficientsSecond),1);
    if(isBasePointIncluded)
        currentFunctionValues(1) = functionValues(1);
        perturbationIndex = 2;
    else
        perturbationIndex = 1;
    end

    for iVariable = 1:nDimension
        rangeIndex = [currentIndex:(currentIndex + nPerturbationSecond - 1)]';

        % Extract function values for this variable
        currentFunctionValues(perturbationIndex:end) = functionValues(rangeIndex,:);
        hessianArray(iVariable, iVariable) = (currentFunctionValues' * activeCoefficientsSecond) / (stepSize(iVariable)^2);

        currentIndex = currentIndex + nPerturbationSecond;
    end
    
    % Compute off-diagonal terms (mixed derivatives)
    nActiveFirst = length(activeCoefficientsFirst);
    for iVariable = 1:nDimension
        for jVariable = (iVariable+1):nDimension
            % Extract function values for this pair of variables
            pairFunctionValues = nan(nActiveFirst, nActiveFirst);

            for iPerturbation = 1:nActiveFirst
                for jPerturbation = 1:nActiveFirst
                    if(~isBasePointIncluded || (iPerturbation ~= 1 && jPerturbation ~= 1))
                        pairFunctionValues(iPerturbation, jPerturbation) = functionValues(currentIndex);
                        currentIndex = currentIndex + 1;
                    end
                end
            end

            if(isBasePointIncluded)
                iPertubationIndex = 1 + (iVariable-1)*nPerturbationSecond + (1:nPerturbationSecond);
                jPertubationIndex = 1 + (jVariable-1)*nPerturbationSecond + (1:nPerturbationSecond);

                pairFunctionValues(1,2:end) = functionValues(iPertubationIndex');
                pairFunctionValues(2:end,1) = functionValues(jPertubationIndex');

                pairFunctionValues(1,1) = functionValues(1);
            end
            
            % Compute mixed derivative using tensor product of first derivative weights
            mixedDerivative = (activeCoefficientsFirst' * pairFunctionValues * activeCoefficientsFirst) / (stepSize(iVariable)*stepSize(jVariable));
            
            % Hessian is symmetric
            hessianArray(iVariable, jVariable) = mixedDerivative;
            hessianArray(jVariable, iVariable) = mixedDerivative;
        end
    end
end
