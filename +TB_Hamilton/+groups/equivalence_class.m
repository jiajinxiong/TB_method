function [ECs,dECS] = equivalence_class(group)
arguments
    group TB_Hamilton.groups.PointGroupElement;
end
ECs = dictionary;
G0 = group;
GK = group.keyhash; G = group;
while true
    if isempty(G)
        break;
    end
    EC = arrayfun(@(x) x.inv*G(1)*x,G0);
    ECK = EC.keyhash;
    [ECK,id] = unique(ECK); EC = EC(id);
    [~,id,~] = intersect(GK,ECK);
    G(id) = []; GK(id) = [];
    if ~ECs.isConfigured
        ECs = insert(ECs,{EC},1);
    else
        if ~isKey(ECs,{EC})
            ECs = insert(ECs,{EC},1);
        end
    end
end

A = ECs.keys;
Iden_id = find(cellfun(@(x) any(ismember(x,group(1).identity)),A));
A = {A{Iden_id},A{1:Iden_id-1},A{Iden_id+1:end}};
[nECS_,Id] = sort(cellfun(@(x) length(x),A));

ECs = A(Id);
if nargout==2
dECS = nECS_;
end

end