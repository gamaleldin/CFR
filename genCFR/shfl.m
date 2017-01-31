% origTensor: original tensor
% shflAcross: a boolean vector of 3 elements, if first elemeny=t is true
% the shuffle is performed across times, if second element is true the
% shuffle is performed across neurons, if third element is true the shuffle
% is performed across conditions. Multiple shuffles are allowed. for
% example, if second and third elements are both true, the shuffle is
% performed across neurons and conditions.
% cyclicShfl: if true the shuffle is performed cyclicly.
%%              
function shuffledTensor = shfl(origTensor, shfl_mode, fix_mode, cyclicShfl)
shuffledTensor = origTensor;
dims = size(origTensor);
tensorIxs = 1:length(dims);

T = dims(shfl_mode);
N = dims(fix_mode);
reorderIx = [fix_mode, shfl_mode, sort(tensorIxs(tensorIxs ~= fix_mode & tensorIxs ~= shfl_mode))]; 
shuffledTensor = reshape(permute(shuffledTensor, reorderIx), dims(fix_mode), dims(shfl_mode), []);

parfor n = 1:N
    if cyclicShfl
        st = randi(T);
        ts = [st:T, 1:st-1];
    else
        ts = randperm(T);
    end
    A = shuffledTensor(n, :, :);
    shuffledTensor(n, :, :) = A(1, ts, :);
end
shuffledTensor = reshape(shuffledTensor, [N, T, dims(tensorIxs ~= fix_mode & tensorIxs ~= shfl_mode)]);
[~, ix]=sort(reorderIx);
shuffledTensor = permute(shuffledTensor, ix);% put back in the original order

end