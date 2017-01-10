function [f, gradf] = objFnMarginalCovTN(K, varargin)
surrTensor = varargin{1};
targetSigma_T = varargin{2};
targetSigma_N = varargin{3};

[fT, ~, gradfT_unc] = objFnMarginalCovT(K, surrTensor, targetSigma_T);
[fN, ~, gradfN_unc] = objFnMarginalCovN(K, surrTensor, targetSigma_N);
f = (fT+fN)/2;
gradf = projEigSpace((gradfT_unc+gradfN_unc)./2, ones(size(K,2),1));
end
