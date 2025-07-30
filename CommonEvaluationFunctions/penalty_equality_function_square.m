function penaltyEqualityValue = penalty_equality_function_square(equalityConstraintValue,varargin)
%PENALTY_EQUALITY_FUNCTION_SQUARE Compute quadratic penalty for equality constraint violations.
%
%   PENALTYEQUALITYVALUE = PENALTY_EQUALITY_FUNCTION_SQUARE(EQUALITYCONSTRAINTVALUE, ...)
%   evaluates the penalty values for equality constraints by squaring the
%   constraint residuals. This penalizes deviations from zero in the equality
%   constraints.
%
%   INPUTS:
%       EQUALITYCONSTRAINTVALUE    - Vector or matrix of equality constraint residuals
%
%   NAME-VALUE PAIR ARGUMENTS:
%       None currently used, but supported for extensibility.
%
%   OUTPUTS:
%       PENALTYEQUALITYVALUE       - Squared penalty values for equality violations
%                                  (same size as input)
%
%   FUNCTIONALITY:
%       The function applies element-wise squaring to the equality constraint
%       values, thus penalizing any deviation from zero.
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

	penaltyEqualityValue = equalityConstraintValue.^2;
end