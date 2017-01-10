function [f, gradf, gradf_unc] = objFnMarginalCovN(K, varargin)
surrTensor = varargin{1};
targetSigma_N = varargin{2};
normalizeTerm = trace(targetSigma_N'*targetSigma_N);
[T, N, C] = size(surrTensor);

XN = reshape(permute((surrTensor),[1 3 2]), [], N);
XNK = XN*K;
Sigma_N = XNK'*XNK;
ER = (targetSigma_N - Sigma_N);
f = trace(ER*ER')./normalizeTerm;
%% calculate gradient 
XNXNK = (XN'*XN)*K;
gradf_unc = (-4/normalizeTerm)*XNXNK*ER;

gradf = projEigSpace(gradf_unc, ones(size(K,2),1));

end

