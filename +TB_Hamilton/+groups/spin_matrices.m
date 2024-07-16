function S  = spin_matrices(spin)
% SPIN_MATRICES is used to construct the spin matrix
% 
% Parameters
% ------
% spin : spin must be half-integer
% 
% Return
% ------
% S : (2s+1,2s+1,3) matrix
%   sx = S(:,:,1); sy = S(:,:,2); sz = S(:,:,3);
%
% Ref: https://en.wikipedia.org/wiki/Spin_(physics)
arguments
    spin            (1,1)    double;
end
d = 2 * spin + 1;
assert(d-fix(d)==0,'The spin must be half-integer');
S = zeros(d,d,3);
e = sqrt( 2*(spin+1)*(1:d-1)-(1:d-1).*(2:d) )/2;

S(:,:,1) = full(spdiags([[0,e]' [e,0]'],[1,-1],d,d));
S(:,:,2) = full(1i*spdiags([-[0,e]' [e,0]'],[1,-1],d,d));
S(:,:,3) = diag(d-1:-2:-d+1)/2;

end