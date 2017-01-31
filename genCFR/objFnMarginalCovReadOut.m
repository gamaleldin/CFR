function [f, gradf, gradf_unc] = objFnMarginalCovReadOut(K, varargin)
surrTensor = varargin{1};
targetSigma = varargin{2};
normalizeTerm = trace(targetSigma'*targetSigma);
N = size(surrTensor, 1); % read out dimensionality = N^2

X = reshape(surrTensor, N, []).';
XK = X*K;
estSigma = XK'*XK;
ER = (targetSigma - estSigma);
f = trace(ER*ER')./normalizeTerm;
%% calculate gradient 
XXK = (X'*X)*K;
gradf_unc = (-4/normalizeTerm)*XXK*ER;

gradf = projEigSpace(gradf_unc, ones(size(K,2),1));

end

