function [err,sizeOut] = checkInputSizes(varargin)
% CHECKINPUTSIZES Check if the sizes of the input arguments match.
%   [err,sizeOut] = CHECKINPUTSIZES(X,Y,...) checks wheter the sizes of the
%   input arguments X,Y,... match. When an input argument is a scalar, this
%   is ignored.
%
%   Inputs: numerical scalars or arrays

% Author: Simon Geirnaert, KU Leuven, Department of Electrical Engineering
% (ESAT), STADIUS Center for Dynamical Systems, Signal Processing and Data
% Analytics & Department of Neurosciences, ExpORL
% Correspondence: simon.geirnaert@esat.kuleuven.be

% Based on MATLAB-function internal.stats.statsizechkM

err = 0;

% Extract size of each input array or scalar
paramSizes = cellfun(@size,varargin(:),'UniformOutput',false);

% Check matching sizes input per input
sizeOut = [];
for in = 1:length(paramSizes)
    if isequal(paramSizes{in}, [1 1])
        continue;
    end
    if isempty(sizeOut)
        sizeOut = paramSizes{in};
    else
        if ~isequal(paramSizes{in},sizeOut)
            err = 1;
            break;
        end
    end
end

% If all inputs are scalar, put corresponding outputSize
if isempty(sizeOut)
   sizeOut = [1 1]; 
end

end