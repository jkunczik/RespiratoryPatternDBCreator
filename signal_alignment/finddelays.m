function [d, err] = finddelays(sigs)
%FINDDELAYS Find delay between a number of given signals with similar
%content
%
%   Uses Matlabs finddelay function (xcorr) to find the delays between all
%   possible pairs of signals and returns a delay for every signal. The
%   delay is computed as the result of a least squares optimization over
%   all delays.

% for signals with multiple components, only keep the first
if size(sigs{1}, 2) > 1
    sigs = cellfun(@(x) x(:,1), sigs, 'UniformOutput', false);
end

for i=1:length(sigs)
    sigs{i} = fillmissing(sigs{i}, 'linear');
end

% all possible combinations
c = nchoosek(1:length(sigs), 2); 

% coefficient table for all combinations (cols=sigs, rows=combinations)
A = zeros(size(c,1), length(sigs));
A(sub2ind(size(A), 1:size(c,1), c(:,1)')) = 1;
A(sub2ind(size(A), 1:size(c,1), c(:,2)')) = -1;

% find delay for all combinations
b = zeros(1, size(c,1));
for i=1:size(c,1)
    b(i) = finddelay(sigs{c(i, 1)}, sigs{c(i, 2)});
end

% find the overal optimal delay for all signals
[d, ~, err] = lsqr(A, b');
d = round(d - min(d));
end

