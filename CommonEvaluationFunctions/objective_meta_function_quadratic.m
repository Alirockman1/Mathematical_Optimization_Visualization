function objectiveValue = objective_meta_function_quadratic(objectiveValueBase,varargin)
%OBJECTIVE_META_FUNCTION_QUADRATIC Aggregate objective terms quadratically from base objective evaluation.
%
%   OBJECTIVEVALUE = OBJECTIVE_META_FUNCTION_QUADRATIC(OBJECTIVEVALUEBASE, ...)
%   computes a scalar objective value by quadratically aggregating various
%   component terms from the base objective evaluation using fallback logic.
%
%   INPUTS:
%       OBJECTIVEVALUEBASE         - Raw objective function values (matrix: samples × objectives)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       Additional options are passed through to OBJECTIVE_META_FUNCTION_BASE_TERM
%
%   OUTPUT:
%       OBJECTIVEVALUE             - Aggregated scalar objective value (vector: samples × 1)
%
%   FUNCTIONALITY:
%       This meta-objective function combines various interpretations of the
%       objective (e.g., less-or-more, nominal, one-sided) using fallback logic
%       with quadratic aggregation:
%           1. Two-sided objectives (squared sum)
%           2. One-sided upper (squared)
%           3. One-sided lower (squared)
%           4. Nominal targets (squared)
%           5. Less-or-more objectives (squared)
%       The final scalar value is the row-wise sum of the selected squared terms.
%
%   DEPENDENCIES:
%       - objective_meta_function_base_term
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
	
	[objectiveLessOrMore,objectiveNominal,objectiveOneSidedLower,objectiveOneSidedUpper] = ...
		objective_meta_function_base_term(objectiveValueBase,varargin{:});
	objectiveTwoSided = objectiveOneSidedLower.^2 + objectiveOneSidedUpper.^2;

	objectiveTerm = nan(size(objectiveValueBase));
	objectiveTerm(isnan(objectiveTerm)) = objectiveTwoSided(isnan(objectiveTerm));
	objectiveTerm(isnan(objectiveTerm)) = objectiveOneSidedUpper(isnan(objectiveTerm)).^2;
	objectiveTerm(isnan(objectiveTerm)) = objectiveOneSidedLower(isnan(objectiveTerm)).^2;
	objectiveTerm(isnan(objectiveTerm)) = objectiveNominal(isnan(objectiveTerm)).^2;
	objectiveTerm(isnan(objectiveTerm)) = objectiveLessOrMore(isnan(objectiveTerm)).^2;

	objectiveValue = sum(objectiveTerm,2);
end