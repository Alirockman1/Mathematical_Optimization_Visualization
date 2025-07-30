function barrierValue = barrier_function_inverse(inequalityConstraintValue,varargin)
%BARRIER_FUNCTION_INVERSE Compute inverse barrier penalty for inequality constraints.
%
%   BARRIERVALUE = BARRIER_FUNCTION_INVERSE(INEQUALITYCONSTRAINTVALUE, ...)
%   computes the barrier penalty values for inequality constraints using the
%   inverse function applied to the negated constraint values.
%
%   INPUTS:
%       INEQUALITYCONSTRAINTVALUE  - Vector or matrix of inequality constraint values
%                                    (must be strictly negative to make sense)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       None currently used, but supported for extensibility.
%
%   OUTPUTS:
%       BARRIERVALUE               - Barrier penalty values computed as -1/x
%                                    (same size as input)
%
%   FUNCTIONALITY:
%       The function applies the inverse barrier method, which penalizes
%       constraint values approaching zero from the negative side, enforcing strict
%       inequality constraints.
%
%   DEPENDENCIES:
%       None
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

	barrierValue = -1./(inequalityConstraintValue);
end