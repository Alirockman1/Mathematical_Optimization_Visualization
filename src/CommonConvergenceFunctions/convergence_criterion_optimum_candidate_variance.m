function hasConverged = convergence_criterion_optimum_candidate_variance(...
    optimumCandidate,optimumCandidatePrevious,...
    optimumObjectiveValue,optimumObjectiveValuePrevious,...
    population,populationPrevious,...
    objectiveValuePopulation,objectiveValuePopulationPrevious,...
    varargin)

%CONVERGENCE_CRITERION_OPTIMUM_CANDIDATE_VARIANCE Check convergence based on optimum candidate objective value change.
%
%   HASCONVERGED = CONVERGENCE_CRITERION_OPTIMUM_CANDIDATE_VARIANCE(...
%       OPTIMUMCANDIDATE, OPTIMUMCANDIDATEPREVIOUS,...
%       OPTIMUMOBJECTIVEVALUE, OPTIMUMOBJECTIVEVALUEPREVIOUS,...
%       POPULATION, POPULATIONPREVIOUS,...
%       OBJECTIVEVALUEPOPULATION, OBJECTIVEVALUEPOPULATIONPREVIOUS,...
%       'Name', Value, ...)
%
%   evaluates whether the optimization has converged based on the change in
%   optimum objective value between successive iterations.
%
%   INPUTS:
%       OPTIMUMCANDIDATE               - Best design vector in current iteration
%       OPTIMUMCANDIDATEPREVIOUS       - Best design vector in previous iteration
%       OPTIMUMOBJECTIVEVALUE          - Objective value of current best design
%       OPTIMUMOBJECTIVEVALUEPREVIOUS  - Objective value of previous best design
%       POPULATION                     - Current population matrix
%       POPULATIONPREVIOUS             - Previous population matrix
%       OBJECTIVEVALUEPOPULATION       - Objective values of current population
%       OBJECTIVEVALUEPOPULATIONPREVIOUS - Objective values of previous population
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'Tolerance'                    - Convergence threshold (default: 1e-5)
%
%   OUTPUT:
%       HASCONVERGED                   - Logical value indicating whether convergence has occurred
%
%   FUNCTIONALITY:
%       This criterion checks if the absolute difference between the best
%       objective values in two successive iterations falls below a given tolerance.
%       If either the previous optimum candidate or objective value is empty,
%       convergence is considered not yet achieved.
%
%   DEPENDENCIES:
%       - None (self-contained)
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Main Author)
%   Copyright 2025 Ali Abbas Kapadia (Contributor) 
%   SPDX-License-Identifier: Apache-2.0

	p = inputParser;
	addParameter(p, 'Tolerance', 1e-05, @(x)isnumeric(x)&&isscalar(x));
	parse(p, varargin{:});
	options = p.Results;

    if(isempty(optimumCandidatePrevious) || isempty(optimumObjectiveValuePrevious))
        hasConverged = false;
        return;
    end
    
	variance = abs(optimumObjectiveValue - optimumObjectiveValuePrevious);
	%variance = norm(optimumCandidate - optimumCandidatePrevious)

	hasConverged = (variance<=options.Tolerance);
end
