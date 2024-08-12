function Irr = Irr_decompose(GCT,character_rep,Isnum)
arguments
GCT                    table;
character_rep       (1,:) double;
Isnum               (1,1) logical = false;
end

% CT_G = TB_Hamilton.groups.Character_Table(G);
character_rep = reshape(character_rep,1,[]);
name = GCT.Irr_Rep;
T = GCT.Equiv_Class;

Irr = round(real(T(2:end,:) * diag(T(1,:)) * character_rep'/sum(T(1,:))));

if ~Isnum
    Irr = Irr(find(Irr))+" " + string(name(find(Irr)+1));
end

end