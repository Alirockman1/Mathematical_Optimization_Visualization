function penaltyInequalityValue = penalty_inequality_function_square(inequalityConstraintValue,varargin)
%PENALTY_INEQUALITY_FUNCTION_SQUARE Compute quadratic penalty for inequality constraint violations.
%
%   PENALTYINEQUALITYVALUE = PENALTY_INEQUALITY_FUNCTION_SQUARE(INEQUALITYCONSTRAINTVALUE, ...)
%   evaluates the penalty values for inequality constraints by squaring the
%   positive parts of the constraint violations. Negative or zero values imply
%   no violation and yield zero penalty.
%
%   INPUTS:
%       INEQUALITYCONSTRAINTVALUE  - Vector or matrix of inequality constraint values
%
%   NAME-VALUE PAIR ARGUMENTS:
%       None currently used, but supported for extensibility.
%
%   OUTPUTS:
%       PENALTYINEQUALITYVALUE     - Squared penalty values for constraint violations,
%                                    zero for non-violations (same size as input)
%
%   FUNCTIONALITY:
%       The function applies max(0, x)^2 element-wise, ensuring that only
%       positive constraint violations contribute to the penalty.
%
%   DEPENDENCIES:
%       None
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Main Author)
%	Copyright 2025 Ali Abbas Kapadia (Contributor)
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

	penaltyInequalityValue = max(0,inequalityConstraintValue).^2;
end