function objectiveValue = objective_meta_function_worst_case(objectiveValueBase,varargin)
%OBJECTIVE_META_FUNCTION_WORST_CASE Aggregate objective terms using worst case from base objective evaluation.
%
%   OBJECTIVEVALUE = OBJECTIVE_META_FUNCTION_WORST_CASE(OBJECTIVEVALUEBASE, ...)
%   computes a scalar objective value by selecting the worst case (maximum)
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
%       with worst case aggregation:
%           1. Two-sided objectives (maximum of lower and upper)
%           2. One-sided upper
%           3. One-sided lower
%           4. Nominal targets
%           5. Less-or-more objectives
%       The final scalar value is the row-wise maximum of the selected terms.
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
	objectiveTwoSided = max(objectiveOneSidedLower,objectiveOneSidedUpper);

	objectiveTerm = nan(size(objectiveValueBase));
	objectiveTerm(isnan(objectiveTerm)) = objectiveTwoSided(isnan(objectiveTerm));
	objectiveTerm(isnan(objectiveTerm)) = objectiveOneSidedUpper(isnan(objectiveTerm));
	objectiveTerm(isnan(objectiveTerm)) = objectiveOneSidedLower(isnan(objectiveTerm));
	objectiveTerm(isnan(objectiveTerm)) = objectiveNominal(isnan(objectiveTerm));
	objectiveTerm(isnan(objectiveTerm)) = objectiveLessOrMore(isnan(objectiveTerm));

	objectiveValue = max(objectiveTerm,[],2);
end