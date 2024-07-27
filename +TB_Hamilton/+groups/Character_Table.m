function GCT = Character_Table(G)
% CHARACTER_TABLE - Construct the character table of a group.
%
% Input:
% ------
% G - a group object.
% Output:
% -------
% GCT - a table with the character table of the group.
arguments
    G (1,:)  TB_Hamilton.groups.PointGroupElement;
end


GEC = TB_Hamilton.groups.equivalence_class(G);
EqC = GEC.keys; numEqC = GEC.values;
M = M_matrices(EqC);
nE = length(numEqC);

while true
    M1 = reshape(tensorprod(rand(1,nE),M,2,1),nE,nE)';
    % [V,E,W] = eig(M1); W1 = W';
    [V,E] = eig(M1); W1 = V.';
    E = diag(E);
    if length(unique(E))==length(E)
        break;
    end
end
f = @(x) fun(x,W1,numEqC);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');

X = fsolve(f,rand(nE,2),options);
T = diag(X(:,1))*W1*diag(X(:,2));


T1 = diag(exp(-1i*angle(T(:,1))))*T;
id = find(all(abs(real(T1)-1)<1e-5,2));

T1 = round(real(T1),4)+1i*round(imag(T1),4);
T2 = [T1(id,:);sortrows(T1([1:id-1,id+1:end],:),1)];
Head = cellfun(@(x) x(1),EqC);
if ~isempty(Head(1).U)
    % F = arrayfun(@(x) round(x.U(1,1)),Head);
    % Str_head = string(Head.Tag);
    % 
    % Str_head = Str_head(:,1);
    % Str_head(abs(F-1)>1e-4) = "d" + Str_head(abs(F-1)>1e-4);
    GCT = array2table([numEqC(:)';T2]);
else
    Str_head = string(Head.Tag);
    GCT = array2table([numEqC(:)';T2],VariableNames=cellstr(Str_head));
end

GCT.Properties.DimensionNames = {'Irr_Rep','Equiv_Class'};
GCT.Properties.RowNames=["Mul","A"+(1:nE)];
end



function M = M_matrices(Eq_classes)
% M_MATRICES - Construct the M matrices for the character table.
% ref: Eq. (7.1) in "Handbook of Computational Group Theory" by Derek F. Holt, Bettina Eick, Eamonn A. O'Brien.
nEqC = length(Eq_classes);
M = zeros(nEqC,nEqC,nEqC);
Key_Eq = cellfun(@(x) x.keyhash,Eq_classes,'UniformOutput',false);
for j1 = 1:nEqC
    g_ell = Eq_classes{j1}(1);
    for j2 = 1:nEqC
        g_ell_k = arrayfun(@(x) x.inv * g_ell,Eq_classes{j2});
        key_g_ell_k = g_ell_k.keyhash;
        Id = cellfun(@(x) nnz(ismember(x,key_g_ell_k)),Key_Eq);
        % M(:,j1,j2) = Id;
        M(j2,j1,:) = Id;
    end
end
end

function y = fun(x,T,numEqC)
T1 = diag(x(:,1))*T*diag(x(:,2));


numEqC = diag(numEqC);
y = [T1*numEqC*T1'-sum(diag(numEqC))*eye(size(numEqC,1));T1'*T1-sum(diag(numEqC))*inv(numEqC)];
end