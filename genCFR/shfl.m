% origTensor: original tensor
% shflAcross: a boolean vector of 3 elements, if first elemeny=t is true
% the shuffle is performed across times, if second element is true the
% shuffle is performed across neurons, if third element is true the shuffle
% is performed across conditions. Multiple shuffles are allowed. for
% example, if second and third elements are both true, the shuffle is
% performed across neurons and conditions.
% cyclicShfl: if true the shuffle is performed cyclicly.
%%              
function shuffledTensor = shfl(origTensor, shfl_mode, cyclicShfl)
[T, N, C] = size(origTensor);
shuffledTensor = origTensor;

%% shuffle across times
if shfl_mode == 1;
    parfor n = 1:N
        if cyclicShfl
            startTime = randi(T);
            ts = [startTime:T, 1:startTime-1];
        else
            ts = randperm(T);
        end
        A = shuffledTensor(:, n, :);
        A = A(ts, 1, :);  
        shuffledTensor(:, n, :) = A;
    end
end
%%
if shfl_mode==2
    parfor c = 1:C
        if cyclicShfl
            startNeu = randi(N);
            ns = [startNeu:N, 1:startNeu-1];
        else
            ns = randperm(N);
        end
        A = shuffledTensor(:, :, c);
        A = A(:, ns, 1);  
        shuffledTensor(:, :, c) = A;
    end
end
%% shuffle across conditions
if shfl_mode ==3
    parfor n = 1:N
        if cyclicShfl
            startCond = randi(C);
            cs = [startCond:C, 1:startCond-1];
        else
            cs = randperm(C);
        end
        A = shuffledTensor(:, n, :);
        A = A(:, 1, cs);  
        shuffledTensor(:, n, :) = A;
    end
end

end