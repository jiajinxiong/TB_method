function PBF= Project_Irr(G,Id_Irr,Basis_function,sym_var)
% PROJECT_IRR is a function that projects a basis function to the irreducible characters of a group.
% The effect of projection operator on an arbitrary function `Basis_function` is given by: 
% $$  \hat{P} = \ell_n/|G| \sum_{g \in G} \chi_n(g)^* . g Basis_function  $$
% where \chi_n(g) is the character of the irreducible representation the group G, 
% and \ell_n is the dimension of the representation.

GCT = TB_Hamilton.groups.Character_Table(G,false);
EC = TB_Hamilton.groups.equivalence_class(G);
ECs = EC.keys;
nG = length(G); GCT = GCT.Equiv_Class;
Dim_Id_Irr = GCT(Id_Irr+1,1);

PBF = Dim_Id_Irr * Basis_function;
for j1 = 2:length(ECs)
    EC = ECs{j1};
    for j2 = 1:length(EC)
        PBF = PBF + conj(GCT(Id_Irr+1,j1)) *EC(j2).apply(Basis_function,sym_var);
    end
end
PBF = PBF * Dim_Id_Irr/nG;
end