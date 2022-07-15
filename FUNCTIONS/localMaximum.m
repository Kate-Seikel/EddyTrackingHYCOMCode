function f = localMaximum(x,minDist)
% function f = localMaximum(x,minDist)
% This function returns the indexes of local maximum in the data x.
% x can be a vector or a matrix
% minDist is the minimum distance between two peaks (local maximas)
% minDist should be a vector in which each argument corresponds to it's
% relevent dimension
% Example:
% x = randn(100,30,10);
% minDist = [10 3 5];
% peak = (x,minDist);

ind=find(isnan(x));
x(ind)=0;
if nargin < 2
    minDist = size(x)/10;
end

dimX = length ( size(x) );
if length(minDist) ~= dimX
    % In case minimum distance isn't defined for all of x dimensions
    % I use the first value as the default for all of the dimensions
    minDist = minDist( ones(dimX,1) );
end

% validity checks
minDist = ceil(minDist);
minDist = max( [minDist(:)' ; ones(1,length(minDist))] );
minDist = min( [minDist ; size(x)] );

se = ones(minDist);
X = imdilate(x,se);
f = find(x == X);