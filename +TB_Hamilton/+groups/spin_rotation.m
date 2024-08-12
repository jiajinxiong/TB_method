function U = spin_rotation(n,s,inversion)
% SPIN_ROTATION constructs the unitary spin rotation matrix for rotation
% specified by the vector with angular momentum 's',given by U =
% exp(-i n * s).
% 
% Parameters
% ------
% n : (1,:) matrix
%   Rotation vector. Its norm is the rotation angle.
% s : (1,1) double 
%   spin. It must equal half-integer
% inversion : (1,1) logical = false
%  if true, the rotation is a improper rotation=exp(1i*pi*s) . Default is false.
%
% Return
% ------
% U : (2s+1,2s+1) matrix
%   spin matrix. its size is (2s+1,2s+1).
arguments
    n    (1,:)  double;
    s    double;
    inversion (1,1) logical = false;
end
spin = s;
assert(~inversion || isscalar(spin), 'When inversion is true, spin must be a scalar');

if isscalar(s)
    s = TB_Hamilton.groups.spin_matrices(s);
end

U = expm(-1i*tensorprod(s,n,3,2));

if inversion && mod(spin,1)==0
    U = expm(1i*pi*spin) * U;
end

end