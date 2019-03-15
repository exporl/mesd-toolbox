function varargout = processInputs(varargin)
% PROCESSINPUTS Check inputs and process them such that they are all
% columns of the same size.
%   [X,Y,...] = PROCESSINPUTS(X,Y,...) checks whether the input sizes match
%   when they are arrays and matches all their sizes. The outputs are the
%   input arguments, stacked in column vectors. If an input is a scalar, it
%   is repeated throughout the column vector.
%
%   Inputs: numerical scalars or arrays

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

% Check if inputs are numerical and/or scalar
chkNum = cellfun(@isnumeric,varargin);
chkScal = cellfun(@isscalar,varargin);
if any(~chkNum & ~chkScal)
    error('Wrong input types.');
end

% Transform inputs to column vectors
varargout = cellfun(@(x)reshape(x,[],1),varargin,'UniformOutput',false);

% Check if the input sizes match
[err,sizeOut] = mtt_utils.checkInputSizes(varargout{:});
if err
    error('Input size mismatch.');
end

% Process scalars into vectors
for in = 1:length(varargout)
    if isscalar(varargout{in})
        varargout{in} = varargout{in}.*ones(sizeOut);
    end
end

end