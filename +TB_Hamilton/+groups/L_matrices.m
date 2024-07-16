function L = L_matrices(d)
% L_MATRICES constructs real space rotation generator matrices in d=2 or d=3 dimension
%
% Parameters
% ------
% d : the space dimension
%
% Return
% ------


arguments
    d (1,1) double {mustBeMember(d,[2,3])} = 2;
end
if d==2
    L = zeros(2,2,1);
    L(:,:,1) = 1i*[0,-1;1,0];
else
    L1 = [0 0 0;0 0 -1;0 1 0]; L2 = [0 0 1;0 0 0;-1 0 0]; L3 = [0 -1 0;1 0 0;0 0 0];
    L = 1i*cat(3,L1,L2,L3);
end
end