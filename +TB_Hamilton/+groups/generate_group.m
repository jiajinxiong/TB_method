function group = generate_group(gens)
    arguments
        gens    TB_Hamilton.groups.PointGroupElement;
    end

    % Initialize the group with the generators
    group = gens;
    oldgroup = gens;

    % Use a set to keep track of unique elements
    groupDict = dictionary();
    groupDict(gens) = 1;

    while true
        newgroup = arrayfun(@(y) arrayfun(@(x) y*x, oldgroup), group,UniformOutput=false);
        newgroup = [newgroup{:}];
        Id = isKey(groupDict,newgroup);
        A = newgroup(~Id);
        groupDict(A) = 1;
        if all(Id)
            break;
        end
        [~,id] = unique(A.keyhash);
        oldgroup = A(id);
    end

    % Convert the set back to an array
    group = keys(groupDict);
    % group = [group{:}];

    % Sort the group elements by their keyhash
    group = sort(group);
end