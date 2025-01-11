function PBF= Project_Irr(G,Id_Irr,Basis_function,sym_var,options)
% PROJECT_IRR is a function that projects a basis function to the irreducible characters of a group.
% The effect of projection operator on an arbitrary function `Basis_function` is given by: 
% $$  \hat{P} = \ell_n/|G| \sum_{g \in G} \chi_n(g)^* . g Basis_function  $$
% where \chi_n(g) is the character of the irreducible representation the group G, 
% and \ell_n is the dimension of the representation.

arguments
    G                   TB_Hamilton.groups.PointGroupElement % Group element
    Id_Irr              (1,1) double             % Index of the irreducible representation
    Basis_function      ;                        % Basis function to be projected 
    sym_var             sym;                     % Symbolic variables
    options.method      (1,1) {mustBeMember(options.method,["basis_fun","Hamilton"])} = "Hamilton";
    options.round_off   (1,1) double = 2;        % Round-off precision
end

% Get the character table for the group G
GCT = TB_Hamilton.groups.Character_Table(G,false);

% Get the equivalence classes of the group G
EC = TB_Hamilton.groups.equivalence_class(G);
ECs = EC.keys;          % Get the label of the equivalence classes, e.g., 'A1', 'A2', 'B1', 'B2', etc.
GCT = GCT.Equiv_Class;  % Get the character table for the equivalence classes


% Get the dimension of the irreducible representation
Dim_Id_Irr = GCT(Id_Irr+1,1);


PBF = 0;
for j1 = 1:length(ECs)
    EC = ECs{j1};
    % Loop over each element in the equivalence class
    for j2 = 1:length(EC)
        % Apply the projection operator
        PBF = PBF + conj(GCT(Id_Irr+1,j1)) *EC(j2).apply(Basis_function,"sym_var",sym_var,"method",options.method);
    end
end


% Normalize the projected basis function
nG = length(G); 
PBF = expand(PBF * Dim_Id_Irr/nG);
PBF = vpa(PBF,3);
for j1 = 1:numel(PBF)
    PBF(j1) = mapSymType(PBF(j1), 'constant', @(x) round(x,options.round_off));
end