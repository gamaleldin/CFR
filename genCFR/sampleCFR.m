function [surrTensor, f] = sampleCFR(dataTensor, params)
shflAcross = [0 0 1]; % [times, neurons conditions]
cyclicShfl = true;
surrType = 'TNC';
targetSigmaT = params.margCov{1};
targetSigmaN = params.margCov{2};
targetSigmaC = params.margCov{3};
meanTensor = params.meanTensor;
if exist('params','var')
    if isfield(params, 'surrType')
        surrType = params.surrType;
    end
    if isfield(params, 'cyclicShfl')
        cyclicShfl = params.cyclicShfl;
    end
    if isfield(params, 'shflAcross')
        shflAcross = params.shflAcross;
    end
end
[T, N, C] = size(dataTensor);
%% Generate surrogates
switch surrType
    case 'T'
    %%%%%%% generate surrogate-T
    surrTensor0 = shfl(dataTensor(:,:, randperm(C)), shflAcross, cyclicShfl);
    [surrTensor, f] = optMarginalCovT(surrTensor0, targetSigmaT);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'TN'
    surrTensor0 = shfl(dataTensor(:,:, randperm(C)), shflAcross, cyclicShfl);
    %%%%%%% generate surrogate-TN
    [surrTensor, f] = optMarginalCovTN(surrTensor0, targetSigmaT , targetSigmaN);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'TNC'
    surrTensor0 = shfl(dataTensor, shflAcross, cyclicShfl);
    %%%%%%% generate surrogate-TNC
    [surrTensor, f] = optMarginalCovTNC(surrTensor0, targetSigmaT , targetSigmaN, targetSigmaC);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
surrTensor = surrTensor+meanTensor;
fprintf('cost (initial->final): %.4f -> %.4f \n', f(1), f(end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
