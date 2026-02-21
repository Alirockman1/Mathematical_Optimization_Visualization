function [objectiveValue,evaluationOutput] = evaluate_optimization_objective(objectiveFunction,designSample,varargin)
%EVALUATE_OPTIMIZATION_OBJECTIVE Evaluate objective function with constraints and penalty/barrier methods.
%
%   [OBJECTIVEVALUE, EVALUATIONOUTPUT] = EVALUATE_OPTIMIZATION_OBJECTIVE(OBJECTIVEFUNCTION, DESIGNSAMPLE, ...)
%   evaluates the provided OBJECTIVEFUNCTION at one or more design points DESIGNSAMPLE, applying optional constraint
%   functions, penalty and barrier terms, and a meta-objective function. Returns scalar objective values and optionally
%   a structured evaluation log.
%
%   [...] = EVALUATE_OPTIMIZATION_OBJECTIVE(..., 'Name', Value, ...) accepts additional name-value pair arguments:
%
%       'ObjectiveFunctionOptions'          - Cell array of additional parameters for OBJECTIVEFUNCTION
%       'ObjectiveMetaFunction'             - Function handle for reducing vector objectives (default: @objective_meta_function_linear)
%       'ObjectiveMetaFunctionOptions'      - Cell array of options for the meta function
%       'ConstraintFunctions'               - Function handle that returns both inequality and equality constraints
%       'ConstraintFunctionsOptions'        - Cell array of parameters for ConstraintFunctions
%       'InequalityConstraintFunction'      - Function handle for inequality constraints
%       'InequalityConstraintFunctionOptions' - Cell array of options for inequality constraints
%       'EqualityConstraintFunction'        - Function handle for equality constraints
%       'EqualityConstraintFunctionOptions' - Cell array of options for equality constraints
%       'PenaltyConstant'                   - Penalty coefficient λ (default: 0)
%       'InequalityPenaltyFunction'         - Function handle to penalize inequality violations (default: @penalty_inequality_function_square)
%       'InequalityPenaltyOptions'          - Cell array of options for inequality penalty function
%       'EqualityPenaltyFunction'           - Function handle to penalize equality violations (default: @penalty_equality_function_square)
%       'EqualityPenaltyOptions'            - Cell array of options for equality penalty function
%       'InequalityBarrierFunction'         - Function handle for barrier term on inequality constraints (default: @barrier_function_natural_log)
%       'InequalityBarrierOptions'          - Cell array of options for barrier function
%       'AugmentedLagrangeInequalityMultiplier' - Scalar or vector multiplier for augmented Lagrangian method
%       'AugmentedLagrangeEqualityMultiplier'   - Scalar or vector multiplier for augmented Lagrangian method
%       'EvaluationOutputFormat'            - Output format: 'struct' (default) or 'array'
%
%   OUTPUTS:
%       OBJECTIVEVALUE     - Final scalar objective values (with penalties/barriers applied)
%       EVALUATIONOUTPUT   - Struct or struct array (if 'struct') or flat array (if 'array') with fields:
%           .DesignPoint
%           .ObjectiveValueBase
%           .ObjectiveValueUnconstrained
%           .InequalityConstraintValue
%           .EqualityConstraintValue
%           .PenaltyInequalityValue
%           .BarrierInequalityValue
%           .PenaltyEqualityValue
%
%   FUNCTIONALITY:
%       This function evaluates the raw objective function and adds penalty or barrier terms for any constraint
%       violations. It supports meta-objective reduction, flexible constraint definitions, and optional use of
%       augmented Lagrangian multipliers. The results can be returned in detailed struct format or as flat arrays.
%
%   DEPENDENCIES:
%       - get_array_subset_if_nonempty
%       - objective_meta_function_linear (default)
%       - penalty_inequality_function_square (default)
%       - penalty_equality_function_square (default)
%       - barrier_function_natural_log (default)
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
	parser.addParameter('ObjectiveFunctionOptions',{});
	parser.addParameter('ObjectiveMetaFunction',@objective_meta_function_linear);
	parser.addParameter('ObjectiveMetaFunctionOptions',{});
	parser.addParameter('ConstraintFunctions',[]);
	parser.addParameter('ConstraintFunctionsOptions',{});
	parser.addParameter('InequalityConstraintFunction',[]);
	parser.addParameter('InequalityConstraintFunctionOptions',{});
	parser.addParameter('EqualityConstraintFunction',[]);
	parser.addParameter('EqualityConstraintFunctionOptions',{});
	parser.addParameter('PenaltyConstant',0);
	parser.addParameter('InequalityPenaltyFunction',@penalty_inequality_function_square);
	parser.addParameter('InequalityPenaltyOptions',{});
	parser.addParameter('EqualityPenaltyFunction',@penalty_equality_function_square);
	parser.addParameter('EqualityPenaltyOptions',{});
	parser.addParameter('InequalityBarrierFunction',@barrier_function_natural_log);
	parser.addParameter('InequalityBarrierOptions',{});
	parser.addParameter('AugmentedLagrangeInequalityMultiplier',0);
	parser.addParameter('AugmentedLagrangeEqualityMultiplier',0);
	parser.addParameter('EvaluationOutputFormat','struct');
	parser.parse(varargin{:});
	options = parser.Results;

	% base evaluation
	nSample = size(designSample,1);

	% evaluate objective function
	objectiveValueBase = objectiveFunction(designSample,options.ObjectiveFunctionOptions{:});

	% objective meta function
	if(size(objectiveValueBase,2)>1 && ~isempty(options.ObjectiveMetaFunction))
		objectiveValueUnconstrained = options.ObjectiveMetaFunction(objectiveValueBase,options.ObjectiveMetaFunctionOptions{:});
	else
		objectiveValueUnconstrained = objectiveValueBase;
	end

	% equality and inequality constraints
	inequalityConstraintValue = [];
	equalityConstraintValue = [];
	if(~isempty(options.ConstraintFunctions))
		[inequalityConstraintValue, equalityConstraintValue] = evaluate_optimization_constraint(options.ConstraintFunctions, designSample, options.ConstraintFunctionsOptions{:});
	end
	if(isempty(inequalityConstraintValue) && ~isempty(options.InequalityConstraintFunction))
		inequalityConstraintValue = options.InequalityConstraintFunction(designSample,options.InequalityConstraintFunctionOptions{:});
	end
	if(isempty(equalityConstraintValue) && ~isempty(options.EqualityConstraintFunction))
		equalityConstraintValue = options.EqualityConstraintFunction(designSample,options.EqualityConstraintFunctionOptions{:});
	end

	% penalty/barrier method - inequality
	penaltyInequalityValue = zeros(nSample,1);
	barrierValue = zeros(nSample,1);
	if(~isempty(inequalityConstraintValue))
		if(~isempty(options.InequalityPenaltyFunction))
			penaltyInequalityValue = options.InequalityPenaltyFunction(inequalityConstraintValue,options.InequalityPenaltyOptions{:});
		elseif(~isempty(options.InequalityBarrierFunction))
			barrierValue = options.InequalityBarrierFunction(inequalityConstraintValue,options.InequalityBarrierOptions{:});
		end
	end

	% penalty - equality
	penaltyEqualityValue = zeros(nSample,1);
	if(~isempty(equalityConstraintValue) && ~isempty(options.EqualityPenaltyFunction))
		penaltyEqualityValue = options.EqualityPenaltyFunction(equalityConstraintValue,options.EqualityPenaltyOptions{:});
	end

	% final congregation
	objectiveValue = objectiveValueUnconstrained;

	if(~isempty(inequalityConstraintValue))
		objectiveValue = objectiveValue + ...
			sum(options.AugmentedLagrangeInequalityMultiplier.*inequalityConstraintValue,2);
	end
	if(~isempty(equalityConstraintValue))
		objectiveValue = objectiveValue + ...
			sum(options.AugmentedLagrangeEqualityMultiplier.*equalityConstraintValue,2);
	end

	if(options.PenaltyConstant>0) 
		objectiveValue = objectiveValue + ...
			options.PenaltyConstant * sum(barrierValue,2) + ...
			(1/options.PenaltyConstant) * (sum(penaltyInequalityValue,2) + sum(penaltyEqualityValue,2));
	end


	if(nargout>1)
		for iSample = 1:nSample
            evaluationOutput(iSample) = struct(...
                'DesignPoint',designSample(iSample,:),...
                'ObjectiveValueBase',objectiveValueBase(iSample,:),...
                'ObjectiveValueUnconstrained',objectiveValueUnconstrained(iSample),...
                'InequalityConstraintValue',get_array_subset_if_nonempty(inequalityConstraintValue,iSample),...
                'EqualityConstraintValue',get_array_subset_if_nonempty(equalityConstraintValue,iSample),...
                'PenaltyInequalityValue',get_array_subset_if_nonempty(penaltyInequalityValue,iSample),...
                'BarrierInequalityValue',get_array_subset_if_nonempty(barrierValue,iSample),...
                'PenaltyEqualityValue',get_array_subset_if_nonempty(penaltyEqualityValue,iSample));
		end

		if(strcmpi(options.EvaluationOutputFormat,'array'))
			outputFieldNames = fieldnames(evaluationOutput);
			for i=1:length(outputFieldNames)
				arrayOutput.(outputFieldNames{i}) = vertcat(evaluationOutput.(outputFieldNames{i}));
			end
			evaluationOutput = arrayOutput;
		end
	end
end
