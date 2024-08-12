function GEC = equivalence_class(group)
arguments
    group TB_Hamilton.groups.PointGroupElement;
end
ECs = dictionary({group(1).identity},1);
G0 = group;
GK = group.keyhash; G = group;
while ~isempty(G)
    EC = arrayfun(@(x) x.inv*G(1)*x,G0);
    ECK = EC.keyhash;
    [ECK,id] = unique(ECK); EC = EC(id);
    [~,id,~] = intersect(GK,ECK);
    G(id) = []; GK(id) = [];
    if ~isKey(ECs,{EC})
        ECs = insert(ECs,{sort(EC)},length(EC));
    end
end

A = ECs.keys;

[nECS_,Id] = sort(values(ECs));

GEC = dictionary(A(Id),nECS_);

end