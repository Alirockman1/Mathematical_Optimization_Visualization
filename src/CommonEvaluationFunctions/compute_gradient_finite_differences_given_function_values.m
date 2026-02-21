function gradientArray = compute_gradient_finite_differences_given_function_values(functionValues, isBasePointIncluded, activeCoefficients, stepSize)
%COMPUTE_GRADIENT_FINITE_DIFFERENCES_GIVEN_FUNCTION_VALUES Compute gradient from pre-evaluated function values
%   This function computes the gradient using pre-evaluated function values at
%   perturbation points. It is primarily used as a helper function by other
%   optimization routines that need to compute gradients efficiently from
%   existing function evaluations.
%
%   GRADIENTARRAY = COMPUTE_GRADIENT_FINITE_DIFFERENCES_GIVEN_FUNCTION_VALUES(FUNCTIONVALUES,
%   ISBASEPOINTINCLUDED, ACTIVECOEFFICIENTS, STEPSIZE) computes the gradient
%   using pre-evaluated function values and finite difference coefficients.
%
%   Inputs:
%       - FUNCTIONVALUES : (nPoints, nFunctionOutput) double
%           Pre-evaluated function values at perturbation points
%       - ISBASEPOINTINCLUDED : logical
%           Whether the base point (unperturbed) evaluation is included
%       - ACTIVECOEFFICIENTS : (nCoefficients, 1) double
%           Finite difference coefficients for gradient computation
%       - STEPSIZE : (nSample, nDimension) double
%           Step sizes used for perturbations in each dimension
%
%   Output:
%       - GRADIENTARRAY : (nFunctionOutput, nDimension, nSample) double
%           Computed gradient array
%
%   See also compute_gradient_finite_differences, compute_hessian_given_function_values,
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
        gradientArray = [];
        return;
    end

    % Extract and reshape values
    nFunctionOutput = size(functionValues,2);
    nSample = size(stepSize,1);
    nDimension = size(stepSize,2);

    nPerturbation = length(activeCoefficients);
    if(isBasePointIncluded)
        nPerturbation = nPerturbation - 1;
    end

    gradientArray = nan(nFunctionOutput, nDimension, nSample);
    blockSize = nPerturbation*nDimension;
    if(isBasePointIncluded)
        blockSize = blockSize + 1;
    end
    
    currentIndex = 1;
    for iDesign = 1:nSample
        rangeIndex = [currentIndex:(currentIndex + blockSize - 1)]';
        currentFunctionValues = functionValues(rangeIndex,:);

        for iOutput = 1:nFunctionOutput
            gradientArray(iOutput,:,iDesign) = compute_gradient_finite_differences_current_point(currentFunctionValues(:,iOutput), isBasePointIncluded, activeCoefficients, stepSize(iDesign,:));
        end

        currentIndex = currentIndex + blockSize;
    end
end

function gradientArray = compute_gradient_finite_differences_current_point(functionValues, isBasePointIncluded, activeCoefficients, stepSize)
    nDimension = size(stepSize,2);
    gradientArray = nan(1,nDimension);

    if(isBasePointIncluded)
        nPerturbation = length(activeCoefficients) - 1;
        currentIndex = 2;
    else
        nPerturbation = length(activeCoefficients);
        currentIndex = 1;
    end

    currentFunctionValues = nan(length(activeCoefficients),1);
    if(isBasePointIncluded)
        currentFunctionValues(1) = functionValues(1);
        perturbationIndex = 2;
    else
        perturbationIndex = 1;
    end

    for i=1:nDimension
        rangeIndex = [currentIndex:(currentIndex + nPerturbation - 1)]';

        currentFunctionValues(perturbationIndex:end) = functionValues(rangeIndex,:);
        gradientArray(i) = (currentFunctionValues' * activeCoefficients) / stepSize(i);

        currentIndex = currentIndex + nPerturbation;
    end
end