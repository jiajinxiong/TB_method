function g = chiral(realspace_dim,options)
arguments
    realspace_dim   (1,1) double;
    options.U                = [];
end
U = options.U;
R = eye(realspace_dim);
g = TB_Hamilton.groups.PointGroupElement(R,conjugate=false,antisymmetry=true,U=U);
end