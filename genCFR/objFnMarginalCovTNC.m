function [f, gradf] = objFnMarginalCovTNC(K, varargin)
surrTensor = varargin{1};
targetSigma_T = varargin{2};
targetSigma_N = varargin{3};
targetSigma_C = varargin{4};
[fT, ~, gradfT_unc] = objFnMarginalCovOther(K, permute(surrTensor, [2, 1, 3]), targetSigma_T);
[fN, ~, gradfN_unc] = objFnMarginalCovReadOut(K, permute(surrTensor, [2, 1, 3]), targetSigma_N);
[fC, ~, gradfC_unc] = objFnMarginalCovOther(K, permute(surrTensor, [2, 3, 1]), targetSigma_C);
% [fT, ~, gradfT_unc] = objFnMarginalCovT(K, surrTensor, targetSigma_T);
% [fN, ~, gradfN_unc] = objFnMarginalCovN(K, surrTensor, targetSigma_N);
% [fC, ~, gradfC_unc] = objFnMarginalCovC(K, surrTensor, targetSigma_C);
f = (fT+fN+fC)/3;
gradf = projEigSpace((gradfT_unc+gradfN_unc+gradfC_unc)./3, ones(size(K,2),1));
end

