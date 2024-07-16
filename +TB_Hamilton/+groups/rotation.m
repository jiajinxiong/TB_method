function g = rotation(angle,axis,inversion,options)
% ROTATION return a rotation operator.
%
% Parameters
% ------
% angle : scalar
%   Rotation angle
% axis : (1,:) matrix (default none)
%   Rotation axis, optional. If not provided, a 2D rotation is generated
%   around the axis normal to the plane. If a 3D vector is provided, a 3D
%   rotation is generated around this axis. Does not need to be normalized
%   to 1.
% inversion : bool (default false)
%   Whether to generate a rotoinversion. By default a proper rotation is
%   returned. Only valid in 3D.
% U : matrix 
%   The unitary action on the Hilbert space.
% spin : 3-dimension tensor or empty or scalar (optional)
%   Spin representation to unitary action. If scalar is provided, it
%   should be half-integer, see `spin_matrices`. Otherwise a 3D tensor must
%   be provided representing 3 components of the angular momentum operator.
%   The unitary action of rotation operator is `U = exp(-i n⋅s)`. In 2D the
%   z axis is assumed to be rotation axis. 

arguments
    angle               (1,1) ;
    axis                (1,:) double = [];
    inversion           (1,1) logical = false;
    options.U           (:,:) double =[];
    options.spin        double = [];
end
spin = options.spin; U = options.U;
assert(isempty(U)|isempty(spin),'Only one of `U` and `spin` may be provided');

if isempty(axis)
    R = TB_Hamilton.groups.spin_rotation(angle,TB_Hamilton.groups.L_matrices(2));
    if ~isempty(spin)
        U = TB_Hamilton.groups.spin_rotation(angle * [0,0,1],spin);
    end
elseif length(axis) == 3
    n = angle * axis/norm(axis);
    R = TB_Hamilton.groups.spin_rotation(n,TB_Hamilton.groups.L_matrices(3));
    if inversion
        R = -R;
    end
    if ~isempty(spin)
        U = TB_Hamilton.groups.spin_rotation(n,spin);
    end
else
    error('`axia` needs to be `empty` or a 3D vector.');
end

g = TB_Hamilton.groups.PointGroupElement(real(R),"U",U);
end