function g = time_reversal(realspace_dim,options)
arguments
    realspace_dim   (1,1) double;
    options.U               (:,:) double =[];
    options.spin            (:,:) double =[];
end
R = eye(realspace_dim); U = options.U; spin = options.spin;
if ~isempty(U)&&~isempty(spin)
    error('Only one of `U` and `spin` may be provided.')
end
if ~isempty(spin)
    U = TB_Hamilton.groups.spin_rotation([0,pi,0],spin);
end
g = TB_Hamilton.groups.PointGroupElement(R,"conjugate",true,"antisymmetry",false,"U",U);
end