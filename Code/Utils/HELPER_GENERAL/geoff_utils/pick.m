function varargout = pick(r, fh, varargin)
% PICK Picks specific outputs from a return list
% [ret1, ret2, ...] = pick(r, fh, ...)
% 
% r is an array of indices of the desired output argument(s)
% fh is the function handle or function name in a string
% Any additional input arguments are treated as arguments to the called
%     function.
% The desired output argument(s) is (are) returned.
%
% As far as I have been able to determine, Matlab provides no built-in
% functionality for easily getting an output argument that is not the first
% in the output arguments list, without also getting all the preceding
% output arguments. This function is designed to overcome that limitation,
% primarily for times when assigning dummy variables is undesirable or
% impossible (e.g., in anonymous functions).
%
% GKA June 2007

output = cell(max(r),1);
[output{:}] = feval(fh, varargin{:});
varargout = output(r);