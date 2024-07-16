function group = generate_group(gens)

oldgroup = gens;
group = gens;
while true
    newgroup = [];
    for j1 = 1:length(group)
        for j2 = 1:length(oldgroup)
            newgroup = [newgroup,group(j1)*oldgroup(j2)];
        end
    end
    newgroup = setdiff(newgroup,group);
    if ~isempty(newgroup)
        group = [group,newgroup];
        oldgroup = newgroup;
    end
end

end