function [f, gradf, gradf_unc] = objFnMarginalCovT(K, varargin)

surrTensor = varargin{1};
targetSigma_T = varargin{2};
[f, ~, gradf_unc] = objFnMarginalCovOther(K, permute(surrTensor, [2, 1, 3]), targetSigma_T);

% surrTensor = varargin{1};
% targetSigma_T = varargin{2};
% normalizeTerm = trace(targetSigma_T.'*targetSigma_T);
% [T, N, C] = size(surrTensor);
% 
% SN = reshape(permute((surrTensor),[1 3 2]), [], N);
% SNK = SN*K;
% STK = reshape(permute(reshape(SNK', N, T, C), [3, 1, 2]), C*N, T);
% 
% Sigma_T = STK'*STK;
% ER = Sigma_T-targetSigma_T;
% f = trace(ER'*ER)/normalizeTerm;
% %% calculate gradient 
% surrTensorT = surrTensor;
% ERsurrTensorT = reshape(ER*reshape(surrTensorT, T, N*C), T, N, C);
% gradf_unc = 4*reshape(permute(surrTensorT, [1, 3, 2]), T*C, N)'*reshape(permute(ERsurrTensorT, [1, 3, 2]), T*C, N)*K;
% gradf_unc = gradf_unc./normalizeTerm;
gradf = projEigSpace(gradf_unc, ones(size(K,2),1));
end
