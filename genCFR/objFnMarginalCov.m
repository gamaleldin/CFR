%% surrTensor; first dimension is the read out dimension, followed by the
% other optimized dimensions
% The order of tensor modes must match the order of the marginal covs
function [f, gradf] = objFnMarginalCov(K, params)
surrTensor = params.surrTensor;   %
margCov = params.margCov;
normFactor = sum(~cellfun(@isempty, margCov));
f =0;
gradf = 0;
if ~isempty(margCov{1})
   [f, gradf] = objFnReadOut(K, surrTensor, margCov{1});
end
dims = size(surrTensor);

tensorIxs = 2:length(dims); %excluding the readout mode
    
for i = tensorIxs
    if ~isempty(margCov{i}) 
        reorderIx = [1, i, sort(tensorIxs(tensorIxs ~= i))]; % make desidred mode second mode;    
        [fi, gradfi] = objFnOther(K, permute(surrTensor, reorderIx), margCov{i});
        f = f + fi;
        gradf = gradf + gradfi;
    end
end
f = f/normFactor;
gradf = gradf/normFactor;
gradf = projEigSpace(gradf, ones(size(K,2),1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, gradf] = objFnReadOut(K, surrTensor, targetSigma)
normalizeTerm = trace(targetSigma'*targetSigma);
N = size(surrTensor, 1); % read out dimensionality = N^2
X = reshape(surrTensor, N, []).';
XK = X*K;
estSigma = XK'*XK;
ER = (targetSigma - estSigma);
f = trace(ER*ER')./normalizeTerm;
%% calculate gradient 
XXK = (X'*X)*K;
gradf = (-4/normalizeTerm)*XXK*ER;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, gradf] = objFnOther(K, surrTensor, targetSigma)
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
gradf = 4*reshape(permute(surrTensorT, [1, 3, 2]), T*C, N)'*reshape(permute(ERsurrTensorT, [1, 3, 2]), T*C, N)*K;
gradf = gradf./normalizeTerm;

end

