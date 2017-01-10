function [f, gradf, gradf_unc] = objFnMarginalCovC(K, varargin)
surrTensor = varargin{1};
targetSigma_C = varargin{2};
normalizeTerm = trace(targetSigma_C.'*targetSigma_C);
[T, N, C] = size(surrTensor);

SN = reshape(permute((surrTensor),[1 3 2]), [], N);
SNK = SN*K;
SCK = reshape(permute(reshape(SNK', N, T, C), [2 1 3]), N*T, C);
Sigma_C = SCK'*SCK;
ER = Sigma_C-targetSigma_C;
f = trace(ER'*ER)/normalizeTerm;

%% calculate gradient 
surrTensorC = permute(surrTensor, [3 2 1]);
ERsurrTensorC = reshape(ER*reshape(surrTensorC, C, N*T), C, N, T);
gradf_unc = 4*reshape(permute(surrTensorC, [1, 3, 2]), C*T, N)'*reshape(permute(ERsurrTensorC, [1, 3, 2]), C*T, N)*K;
gradf_unc = gradf_unc./normalizeTerm;
gradf = projEigSpace(gradf_unc, ones(size(K,2),1));
end
