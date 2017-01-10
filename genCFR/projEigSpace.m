function [G] = projEigSpace(grad, v)
v = v./norm(v);
N = length(v);
I = eye(N);
G = grad*(I-(v*v'));
end