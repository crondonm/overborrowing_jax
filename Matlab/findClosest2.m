function inds = findClosest2(range, vector)
% ‘flatten’ everything into row vectors
origshape = size(vector);
sz_r = size(range);
numTargets = prod(sz_r);
numInputs = prod(origshape);

range = reshape(range, [1, numTargets]);
vector = reshape(vector, [1, numInputs]);

% find the difference between the target values and the actual values
% ‘range’ is replicated across columns,
% ‘vector’ is replicated across rows
diff2 = repmat(range', 1, numInputs) - repmat(vector, numTargets, 1);
diff4=abs(diff2)';
[~, inds] = min(diff4);
inds = reshape(inds, sz_r);
end