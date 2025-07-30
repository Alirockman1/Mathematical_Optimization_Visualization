function barrierValue = barrier_function_natural_log(inequalityConstraintValue,varargin)
%BARRIER_FUNCTION_NATURAL_LOG Compute natural logarithm barrier penalty for inequality constraints.
%
%   BARRIERVALUE = BARRIER_FUNCTION_NATURAL_LOG(INEQUALITYCONSTRAINTVALUE, ...)
%   computes the barrier penalty values for inequality constraints using the
%   natural logarithm function applied to the negated constraint values.
%
%   INPUTS:
%       INEQUALITYCONSTRAINTVALUE  - Vector or matrix of inequality constraint values
%                                    (must be strictly negative to avoid complex results)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       None currently used, but supported for extensibility.
%
%   OUTPUTS:
%       BARRIERVALUE               - Barrier penalty values computed as -log(-x)
%                                    (same size as input)
%
%   FUNCTIONALITY:
%       The function applies the natural logarithm barrier method, which penalizes
%       constraint values approaching zero from the negative side, enforcing strict
%       inequality constraints.
%
%   NOTES:
%       The input values must be strictly less than zero; otherwise, the logarithm
%       is undefined or complex.
%
%   DEPENDENCIES:
%       None
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Main Author)
%   Copyright 2025 Ali Abbas Kapadia (Contributor)
%   SPDX-License-Identifier: Apache-2.0

	barrierValue = -log(-inequalityConstraintValue);
end