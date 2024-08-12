function GCT = Character_Table(G,display)
% CHARACTER_TABLE - Construct the character table of a group.
% 1D-A,B; 2D-E; 3D-T,F; 4D-G;
%
% Input:
% ------
% G - a group object.
% Output:
% -------
% GCT - a table with the character table of the group.
arguments
    G TB_Hamilton.groups.PointGroupElement;
    display (1,1) logical = true;
end


GEC = TB_Hamilton.groups.equivalence_class(G);
EqC = GEC.keys; numEqClasses = GEC.values;
M = M_matrices(EqC);
nE = length(numEqClasses);

while true
    M1 = reshape(tensorprod(rand(1,nE),M,2,1),nE,nE)';
    % [V,E,W] = eig(M1); W1 = W';
    [V,E] = eig(M1); W1 = V.';
    E = diag(E);
    if length(unique(E))==length(E)
        break;
    end
end



f = @(x) fun(x,W1,numEqClasses);
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','none');

% X = fsolve(f,rand(nE,2),options);
% T = diag(X(:,1))*W1*diag(X(:,2));

X = fsolve(f,rand(nE,1),options);
T = diag(X)*W1*inv(diag(numEqClasses));


T1 = diag(exp(-1i*angle(T(:,1))))*T;
id = find(all(abs(real(T1)-1)<1e-5,2));

T1 = round(real(T1),4)+1i*round(imag(T1),4);
T2 = [T1(id,:);sortrows(T1([1:id-1,id+1:end],:))];
Head = cellfun(@(x) x(1),EqC);   % A representation of the equivalence class.


head_name = equiv_class_name(Head);
GCT = array2table([numEqClasses(:)';T2],VariableNames=cellstr(head_name));
GCT.Properties.DimensionNames = {'Irr_Rep','Equiv_Class'};

GCT.Properties.RowNames=Irr_name(T2(:,1));

if display
    fig = uifigure("Position",[500 500 760 360]);
    uitTable = uitable(fig, ...
        "Data",GCT, ...
        "Position",[20 20 720 320]);
    % % Center all cells of the table & bolden the first row
    uisCenter = uistyle('HorizontalAlignment', 'center','FontName','Arial');
    % uisBold = uistyle('FontWeight','bold');
    addStyle(uitTable, uisCenter)
end


end

function y = fun(x,T,numEqC)
% T1 = diag(x(:,1))*T*diag(x(:,2));
T1 = diag(x)*T*inv(diag(numEqC));
numEqC1 = diag(numEqC);
y = [T1*numEqC1*T1'-sum(diag(numEqC1))*eye(size(numEqC1,1));T1'*T1-sum(diag(numEqC1))*inv(numEqC1)]; %#ok<MINV>
end

function M = M_matrices(Eq_classes)
% M_MATRICES - Construct the M matrices for the character table.
% ref: Eq. (7.1) in "Handbook of Computational Group Theory" by Derek F. Holt, Bettina Eick, Eamonn A. O'Brien.
%
% Input:
% ------
% Eq_classes - a cell array of equivalence classes.
%
% Output:
% -------
% M - a 3D array with the M matrices.

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

function Str_head = equiv_class_name(Head)
arguments
    Head (1,:) TB_Hamilton.groups.PointGroupElement
end
Str_head = string(Head.Tag);
Str_head = Str_head(:,1);
unique_str = unique(Str_head);
if length(Str_head)/length(unique_str)==2
    for j1 = 1:length(unique_str)
        nid = find(Str_head==unique_str(j1));
        Str_head(nid) = ["","d"]'  + Str_head(nid);
    end
else
    for j1 = 1:length(unique_str)
        nid = find(Str_head==unique_str(j1));
        Str_head(nid) = ["","d_"+string(1:length(nid)-1)]'+" "+ Str_head(nid);
    end
end

end

function name = Irr_name(Irr_D)
Irr_D = [0,Irr_D'];
name = strings(1,length(Irr_D));
name(1) = "Mult.";
name(Irr_D==1) = "A"+(1:nnz(Irr_D==1)); name(Irr_D==2) = "E"+(1:nnz(Irr_D==2));
name(Irr_D==3) = "F"+(1:nnz(Irr_D==3)); name(Irr_D==4) = "G"+(1:nnz(Irr_D==4));
name(Irr_D==5) = "H"+(1:nnz(Irr_D==5)); name(Irr_D>5) =  "I"+(1:nnz(Irr_D>5));
end