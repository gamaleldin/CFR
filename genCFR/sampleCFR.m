function [surrTensor, f] = sampleCFR(dataTensor, params)
readout_mode = 1;
shfl_mode = 3;
fix_mode = 2;
cyclicShfl = true;
margCov = params.margCov;
meanTensor = params.meanTensor;
dims = size(dataTensor);
tensorIxs = 1:length(dims);
if exist('params','var')
    if isfield(params, 'cyclicShfl')
        cyclicShfl = params.cyclicShfl;
    end
    if isfield(params, 'shfl_mode')
        shfl_mode = params.shfl_mode;
    end
    if isfield(params, 'readout_mode')
        readout_mode = params.readout_mode;
    end
    
    if isfield(params, 'fix_mode')
        fix_mode = params.fix_mode;
    end
end

%% Generate surrogates
%% 1- shuffling step
surrTensor0 = shfl(dataTensor, shfl_mode, fix_mode, cyclicShfl);  % shuffle data
%% 2- correction step
% make readout mode first
reorderIx = [readout_mode, sort(tensorIxs(tensorIxs ~= readout_mode))]; 
[surrTensor, f, K] = optMarginalCov(permute(surrTensor0, reorderIx), margCov(reorderIx));
[~, ix]=sort(reorderIx);
surrTensor = permute(surrTensor, ix);% put back in the original order

%% add mean tensor
surrTensor = surrTensor+meanTensor;  
fprintf('cost (initial->final): %.4f -> %.4f \n', f(1), f(end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
