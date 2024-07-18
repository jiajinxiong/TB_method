function group = generate_group(gens)

oldgroup = gens;
group = gens;

while true
    newgroup = [];
    newgroup = arrayfun(@(y) cat(2,arrayfun(@(x) y*x, oldgroup),newgroup), group,UniformOutput=false);
    newgroup = [newgroup{:}];

    % newgroup_tag = arrayfun(@(x) x.tag,newgroup);
    NGT = TB_Hamilton.groups.grouptable(newgroup);
    GT = TB_Hamilton.groups.grouptable(group);
    % group_tag = arrayfun(@(x) x.tag,group);
    % [~,id] = setdiff(newgroup_tag,group_tag);
    [~,id] = setdiff(NGT,GT,"rows");
    if isempty(id)
        break;
    end
    group = cat(2,group,newgroup(id'));
end
% G1 = []; G1_tag = [];
% for j1 = 1:length(group)
%     id = find(group==group(j1),1);
%     G1 = [G1,group(id)];
%     G1_tag = [G1_tag,group_tag(id)];
% end


[~,id] = sortrows(GT);
group = group(id');
end


