function [surrTensorOut, f, K] = optMarginalCov(surrTensor, margCov)
dims = size(surrTensor);
N = dims(1); % readout mode dimensionality
meanN = sumTensor(surrTensor, 2:length(dims))/prod(dims(2:end));
surrTensor = bsxfun(@minus, surrTensor, meanN);
%% solve for K
K = eye(N);
K = K*(eye(N)-ones(N)./N);
maxIter = 100;
params.surrTensor = surrTensor;   %
params.margCov = margCov;

[K, f] = minimize(K ,'objFnMarginalCov' , maxIter, params);

%%
surrTensorOut = reshape(K.'*reshape(surrTensor, N, []), [N, dims(2:end)]);
end

