function [surrTensorOut, f, K] = optMarginalCovTNC(surrTensor, targetSigma_T, targetSigma_N, targetSigma_C)
[T, N, C] = size(surrTensor);
meanN = sumTensor(surrTensor, [1 3])/(C*T);
surrTensor = bsxfun(@minus, surrTensor, meanN);
%% solve for K
K = eye(N);
K = K*(eye(N)-ones(N)./N);
maxIter = 100;
[K, f] = minimize(K ,'objFnMarginalCovTNC' , maxIter, surrTensor , targetSigma_T, targetSigma_N, targetSigma_C);

%%
surrTensorOut = permute(reshape(...
    (reshape(permute(surrTensor,[1 3 2]), [], N)*K).',...
    N, T, C), [2 1 3]);
end



