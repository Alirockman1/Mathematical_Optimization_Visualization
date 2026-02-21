function [inequalityConstraintValue, equalityConstraintValue, isSingleOutput] = evaluate_optimization_constraint(constraintFunction, designSample, varargin)
%EVALUATE_OPTIMIZATION_CONSTRAINT Evaluate constraint function(s) for optimization
%   This function evaluates a constraint function handle at given design points.
%   It supports constraint functions that return either only inequality constraints
%   or both inequality and equality constraints.
%
%   [INEQUALITYCONSTRAINTVALUE, EQUALITYCONSTRAINTVALUE, ISSINGLEOUTPUT] = EVALUATE_OPTIMIZATION_CONSTRAINT(CONSTRAINTFUNCTION, DESIGNSAMPLE)
%   evaluates the constraint function at the provided design points.
%
%   Output behavior:
%     - If the constraint function returns a single output, INEQUALITYCONSTRAINTVALUE contains the result and EQUALITYCONSTRAINTVALUE is empty.
%     - If the constraint function returns two outputs, INEQUALITYCONSTRAINTVALUE and EQUALITYCONSTRAINTVALUE contain the respective results.
%     - ISSINGLEOUTPUT is true if only one output is returned, false otherwise.
%
%   EVALUATE_OPTIMIZATION_CONSTRAINT(..., NAME, VALUE, ...) allows specification of
%   additional options through name-value pairs, which are forwarded to the constraint function.
%
%   Inputs:
%       - CONSTRAINTFUNCTION : function handle
%           Constraint function to evaluate. Must return either one (inequality) or two (inequality, equality) outputs.
%       - DESIGNSAMPLE : (nSample, nDimension) double
%           Design points at which to evaluate the constraint(s).
%       - ... : Additional name-value pairs passed to the constraint function.
%
%   Outputs:
%       - INEQUALITYCONSTRAINTVALUE : array
%           Values of the inequality constraints at the design points.
%       - EQUALITYCONSTRAINTVALUE : array
%           Values of the equality constraints at the design points (empty if not present).
%       - ISSINGLEOUTPUT : logical
%           True if only one output is returned by the constraint function, false otherwise.
%
%   See also evaluate_optimization_objective, compute_gradient_finite_differences, compute_hessian_finite_differences.
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

    isSingleOutput = true;
    equalityConstraintValue = [];

    if(nargout(constraintFunction) == 1)
        inequalityConstraintValue = constraintFunction(designSample,varargin{:});
    elseif(nargout(constraintFunction) == 2)
        [inequalityConstraintValue,equalityConstraintValue] = constraintFunction(designSample,varargin{:});
        isSingleOutput = false;
    else
        try
            [inequalityConstraintValue,equalityConstraintValue] = constraintFunction(designSample,varargin{:});
            isSingleOutput = false;
        catch
            try
                inequalityConstraintValue = constraintFunction(designSample,varargin{:});
            catch
                error('evaluate_optimization_constraint:InvalidNumberOfOutputs', 'The function must return either one or two outputs.');
            end
        end
    end
end