function group = generate_group(gens)

oldgroup = gens;
group = gens;

while true
    newgroup = [];
    for j1 = 1:length(group)
        for j2 = 1:length(oldgroup)
            g0 = group(j1)*oldgroup(j2);
            if ~any([group,newgroup] == g0)
                newgroup = [g0,newgroup];
            end
        end
    end

    if ~isempty(newgroup)
        group = [group,newgroup];
        oldgroup = newgroup;
    else
        break;
    end
end
group = num2cell(group);
end