function maxidx = argmax(input, dim)
%argmax   Indices of largest elements in array

% Check for proper input arguments
if nargin < 1 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

if nargin < 2 || isempty(dim)
    [temp, maxidx] = max(input);
else
    [temp, maxidx] = max(input,[],dim);
end
