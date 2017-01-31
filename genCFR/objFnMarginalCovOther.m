% marginal cov placed in second mode, read out in the first mode of the
% tensor
function [f, gradf, gradf_unc] = objFnMarginalCovOther(K, varargin)
surrTensor = varargin{1};
targetSigma = varargin{2};
normalizeTerm = trace(targetSigma.'*targetSigma);
N = size(surrTensor, 1); % read out size
T = size(surrTensor, 2); % mode of interest size
C = round(numel(surrTensor)/(T*N)); % the combined dimensionality of all other modes
surrTensor =  reshape(surrTensor, N, T, C); % reduce the tensor to 3 mode tensor with first mode is the read out second mode is the desired mode

X = reshape(surrTensor, N, T*C).';
XK = X*K;

XTK = reshape(permute(reshape(XK.', N, T, C), [3, 1, 2]), C*N, T);

estSigma = XTK'*XTK;
ER = estSigma-targetSigma;
f = trace(ER'*ER)/normalizeTerm;
%% calculate gradient 
surrTensorT = permute(surrTensor, [2 1 3]); % make the desired mode first
ERsurrTensorT = reshape(ER*reshape(surrTensorT, T, N*C), T, N, C);
gradf_unc = 4*reshape(permute(surrTensorT, [1, 3, 2]), T*C, N)'*reshape(permute(ERsurrTensorT, [1, 3, 2]), T*C, N)*K;
gradf_unc = gradf_unc./normalizeTerm;
gradf = projEigSpace(gradf_unc, ones(size(K,2),1));

end
