function SG = find_symmetries(H,group,sym_var)
arguments
    H           sym ;
    group       TB_Hamilton.groups.PointGroupElement;
    sym_var     (1,:)   sym;
end

% H = simplify(H,'Steps',200);
Id = []; 
for j1 = 1:length(group)
    Delta = group(j1).apply(H,"sym_var",sym_var)-H;
    % Delta = simplify(Delta);
    sym_vars = symvar(Delta);
    Delta = double(subs(Delta, sym_vars ,rand(1,length(sym_vars))));
    if max(abs(Delta))<1e-6
        Id = [Id,j1];
    end
end
SG = group(Id);
end