function U = spin_rotation(n,s)
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
%
% Return
% ------
% U : (2s+1,2s+1) matrix
%   spin matrix. its size is (2s+1,2s+1).
arguments
    n    (1,:)  double;
    s    double;
end

if length(s) == 1
    s = TB_Hamilton.groups.spin_matrices(s);
end

U = expm(-1i*tensorprod(s,n,3,2));

end