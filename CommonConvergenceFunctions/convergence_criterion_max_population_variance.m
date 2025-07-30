function hasConverged = convergence_criterion_max_population_variance(...
    optimumCandidate,optimumCandidatePrevious,...
    optimumObjectiveValue,optimumObjectiveValuePrevious,...
    population,populationPrevious,...
    objectiveValuePopulation,objectiveValuePopulationPrevious,...
    varargin)
%CONVERGENCE_CRITERION_MAX_POPULATION_VARIANCE Optimization end-condition
%	CONVERGENCE_CRITERION_MAX_POPULATION_VARIANCE is a function that checks if the
%	maximum population variance is less than a given tolerance. This can be used
%   as a convergence criterion for optimization algorithms that use a population
%   of designs.
%
%   HASCONVERGED = CONVERGENCE_CRITERION_MAX_POPULATION_VARIANCE(...
%       OPTIMUMCANDIDATE, OPTIMUMCANDIDATEPREVIOUS,...
%       OPTIMUMOBJECTIVEVALUE, OPTIMUMOBJECTIVEVALUEPREVIOUS,...
%       POPULATION, POPULATIONPREVIOUS,...
%       OBJECTIVEVALUEPOPULATION, OBJECTIVEVALUEPOPULATIONPREVIOUS,...
%       'Name', Value, ...)
%
%   Inputs:
%       - optimumCandidate : (1,nDimension) double
%       - optimumCandidatePrevious : (1,nDimension) double
%       - optimumObjectiveValue : double
%       - optimumObjectiveValuePrevious : double
%       - population : (nPopulation,nDimension) double
%       - populationPrevious : (nPopulation,nDimension) double
%       - objectiveValuePopulation : (nPopulation,1) double
%       - objectiveValuePopulationPrevious : (nPopulation,1) double
%
%   NAME-VALUE PAIR ARGUMENTS:
%       Tolerance : (1,1) double
%           Tolerance for the convergence criterion.
%
%   OUTPUTS:
%       hasConverged : (1,1) logical
%
%   See also convergence_criterion_optimum_candidate_variance.
%	
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Gaurav Vaibhav (Main Author)
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
	addParameter(p, 'Tolerance', 1e-05, @(x)isnumeric(x)&&isscalar(x));
	parse(p, varargin{:});
	options = p.Results;

	maxPopulationValue = max(objectiveValuePopulation);
	minPopulationValue = min(objectiveValuePopulation);
	variance = abs(maxPopulationValue - minPopulationValue);

	hasConverged = (variance<=options.Tolerance);
end
