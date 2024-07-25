function GT = grouptable(G)
arguments
    G TB_Hamilton.groups.PointGroupElement;
end

GT = table;
GT.rotation = arrayfun(@(x) mat2str(round(x.R,2)),G,'UniformOutput',false)';

c = [G.conjugate];
if any(c)
    GT.conjugate = c';
end
a = [G.antisymmetry];
if any(a)
    GT.antisymmetry = a';
end

if ~isempty(G(1).U)
    GT.U = arrayfun(@(x) mat2str(determine_phase(x.U)),G,'UniformOutput',false)';
end

end


function y = determine_phase(x)
x(abs(x)<1e-3) = 0;
id = find(x,1);
y = x * exp(-1i*angle(x(id)/abs(x(id))));
y = round(real(y),2) + round(imag(y),2)*1i;
end