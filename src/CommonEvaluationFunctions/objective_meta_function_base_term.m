function [objectiveLessOrMore,objectiveNominal,objectiveOneSidedLower,objectiveOneSidedUpper] = objective_meta_function_base_term(objectiveValueBase,varargin)
%OBJECTIVE_META_FUNCTION_BASE_TERM Compute base objective terms with normalization and weighting.
%
%   [OBJECTIVELESSORMORE, OBJECTIVENOMINAL, OBJECTIVEONESIDEDLOWER, OBJECTIVEONESIDEDUPPER] = ...
%       OBJECTIVE_META_FUNCTION_BASE_TERM(OBJECTIVEVALUEBASE, ...)
%   computes several forms of objective function components derived from the
%   raw objective values. These terms facilitate flexible objective formulations
%   such as nominal, one-sided, and weighted objectives with normalization.
%
%   INPUTS:
%       OBJECTIVEVALUEBASE         - Raw objective function values (matrix: samples × objectives)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'Weights'                 - Vector of weights for each objective (default: ones)
%       'NormalizationFactor'     - Vector of normalization factors (default: computed from bounds)
%       'NominalValue'            - Nominal target values for objectives (default: NaN)
%       'LowerLimit'              - Lower limits for one-sided objectives (default: NaN)
%       'UpperLimit'              - Upper limits for one-sided objectives (default: NaN)
%
%   OUTPUTS:
%       OBJECTIVELESSORMORE        - Weighted raw objective values (samples × objectives)
%       OBJECTIVENOMINAL           - Weighted absolute deviation from nominal values
%       OBJECTIVEONESIDEDLOWER     - Weighted normalized penalties below lower limit
%       OBJECTIVEONESIDEDUPPER     - Weighted normalized penalties above upper limit
%
%   FUNCTIONALITY:
%       This function calculates four types of objective components:
%       1. Less-or-more objectives (raw weighted values)
%       2. Nominal objectives (absolute deviation from nominal)
%       3. One-sided lower penalties normalized by a factor
%       4. One-sided upper penalties normalized by a factor
%       Normalization factors are determined automatically if not provided,
%       using a hierarchy of available bounds and nominal values.
%
%   DEPENDENCIES:
%       None
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

	parser = inputParser;
	parser.addParameter('Weights',[]);
	parser.addParameter('NormalizationFactor',[]);
	parser.addParameter('NominalValue',[]);
	parser.addParameter('LowerLimit',[]);
	parser.addParameter('UpperLimit',[]);
	parser.parse(varargin{:});
    options = parser.Results;

	nObjective = size(objectiveValueBase,2);

	nominalValue = options.NominalValue;
	if(isempty(nominalValue))
		nominalValue = nan(1,nObjective);
	end

	lowerLimit = options.LowerLimit;
	if(isempty(lowerLimit))
		lowerLimit = nan(1,nObjective);
	end

	upperLimit = options.UpperLimit;
	if(isempty(upperLimit))
		upperLimit = nan(1,nObjective);
	end

	normalizationFactor = options.NormalizationFactor;
	if(isempty(normalizationFactor))
		normalizationFactor = nan(1,nObjective);
	end
	unsuitableFactor = @(entry) [isinf(entry) | isnan(entry) | (entry==0)];
	i = 1;
	notDefined = unsuitableFactor(normalizationFactor);
	while(any(notDefined))
		switch i
			case 1
				normalizationFactor(notDefined) = upperLimit(notDefined) - lowerLimit(notDefined);
			case 2
				normalizationFactor(notDefined) = abs(upperLimit(notDefined));
			case 3
				normalizationFactor(notDefined) = abs(lowerLimit(notDefined));
			case 4
				normalizationFactor(notDefined) = abs(nominalValue(notDefined));
			otherwise
				normalizationFactor(notDefined) = 1;
		end
		notDefined = unsuitableFactor(normalizationFactor);
		i = i+1;
	end

	weight = options.Weights;
	if(isempty(weight))
		weight = ones(1,nObjective);
	end
	
	objectiveLessOrMore = weight.*objectiveValueBase;
	objectiveNominal = weight.*abs(objectiveValueBase - nominalValue);
	objectiveOneSidedLower = weight.*(lowerLimit - objectiveValueBase)./normalizationFactor;
	objectiveOneSidedUpper = weight.*(objectiveValueBase - upperLimit)./normalizationFactor;
end