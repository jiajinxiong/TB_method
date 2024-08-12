function PBF= Project_Irr(G,Id_Irr,Basis_function,sym_var,options)
% PROJECT_IRR is a function that projects a basis function to the irreducible characters of a group.
% The effect of projection operator on an arbitrary function `Basis_function` is given by: 
% $$  \hat{P} = \ell_n/|G| \sum_{g \in G} \chi_n(g)^* . g Basis_function  $$
% where \chi_n(g) is the character of the irreducible representation the group G, 
% and \ell_n is the dimension of the representation.

arguments
    G                   TB_Hamilton.groups.PointGroupElement
    Id_Irr              (1,1) double
    Basis_function      ;
    sym_var             sym;
    options.method      (1,1) {mustBeMember(options.method,["basis_fun","Hamilton"])} = "Hamilton";
    options.round_off   (1,1) double = 2;
end

GCT = TB_Hamilton.groups.Character_Table(G,false);
EC = TB_Hamilton.groups.equivalence_class(G);
ECs = EC.keys;
GCT = GCT.Equiv_Class;
Dim_Id_Irr = GCT(Id_Irr+1,1);


PBF = 0;
for j1 = 1:length(ECs)
    EC = ECs{j1};
    for j2 = 1:length(EC)
        PBF = PBF + conj(GCT(Id_Irr+1,j1)) *EC(j2).apply(Basis_function,"sym_var",sym_var,"method",options.method);
    end
end
nG = length(G); 
PBF = expand(PBF * Dim_Id_Irr/nG);
PBF = vpa(PBF,3);
for j1 = 1:numel(PBF)
    PBF(j1) = mapSymType(PBF(j1), 'constant', @(x) round(x,options.round_off));
end