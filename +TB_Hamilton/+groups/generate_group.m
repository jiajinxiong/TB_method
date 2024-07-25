function group = generate_group(gens)
arguments
    gens TB_Hamilton.groups.PointGroupElement;
end
[~,id] = unique(gens.keyhash);
gens = gens(id);
oldgroup = gens; group = gens;
GK = group.keyhash;
while true
    newgroup = arrayfun(@(y) arrayfun(@(x) y*x, oldgroup), group,UniformOutput=false);
    newgroup = [newgroup{:}];
    NGK = newgroup.keyhash;
    % GK = group.keyhash;
    [~,id] = setdiff(NGK,GK);
    if isempty(id)
        break;
    end
    group = cat(2,group,newgroup(id'));
    GK = cat(2,GK,NGK(id'));
end
[~,id] = sort(GK); group = group(id');
end