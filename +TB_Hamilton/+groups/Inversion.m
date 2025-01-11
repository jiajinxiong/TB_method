function g = Inversion(realspace_dim,options)
    arguments
        realspace_dim           (1,1) double;
        options.U               (:,:) double =[];
        options.spin            (:,:) double =[];
    end
    R = eye(realspace_dim); U = options.U; spin = options.spin;
    if ~isempty(U)&&~isempty(spin)
        error('Only one of `U` and `spin` may be provided.')
    end
    if ~isempty(spin)
        g = TB_Hamilton.groups.rotation(0,[1,0,0],"inversion",true,"spin",spin);
    elseif ~isempty(U)
        g = TB_Hamilton.groups.rotation(0,[1,0,0],"inversion",true,"U",U);
    else
        g = TB_Hamilton.groups.rotation(0,[1,0,0],"inversion",true);
    end

end