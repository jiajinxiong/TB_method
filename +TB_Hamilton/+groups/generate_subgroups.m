function subgroups = generate_subgroups(group)
arguments
    group TB_Hamilton.groups.PointGroupElement;
end
import TB_Hamilton.groups.generate_group;


g1 = arrayfun(@(x) generate_group(x),group,'UniformOutput',false);
sg1 = dictionary(g1,num2cell(group,length(g1)));

subgroups = sg1; sgold = sg1;
nsg1 = numEntries(sg1);key_sg1 = sg1.keys; value_sg1 = sg1.values;
while true
    sgnew = dictionary;
    nsgold = numEntries(sgold);
    
    for j1 = 1:nsg1
        GK1 = key_sg1{j1}.keyhash;
        for j2 = 1:nsgold
            key_sgold = sgold.keys; value_sgold = sgold.values;
            GK2 = key_sgold{j2}.keyhash;
            if ~isempty(setdiff(GK1,GK2))
                newgen = [value_sg1{j1},value_sgold{j2}];
                newsg = generate_group(newgen);
                if ~isKey(subgroups,{newsg})
                    subgroups = subgroups.insert({newsg},{newgen});
                    sgnew = sgnew.insert({newsg},{newgen});
                end
            end
        end
    end
    if sgnew.numEntries > 0
        sgold = sgnew;
    end
    if sgnew.numEntries==0 || min(cellfun(@(x) length(x),sgnew.keys))==length(group)
        break;
    end
end

subgroups = table(subgroups.keys,subgroups.values,'VariableNames',{'Subgroups','Generator'});

end

% 
% function process_subgroup(j1, j2, key_sg1, value_sg1, key_sgold, value_sgold, subgroups, sgnew)
%     GK1 = key_sg1{j1}.keyhash;
%     GK2 = key_sgold{j2}.keyhash;
%     if ~isempty(setdiff(GK1, GK2))
%         newgen = [value_sg1{j1}, value_sgold{j2}];
%         newsg = generate_group(newgen);
%         if ~isKey(subgroups, {newsg})
%             subgroups.insert({newsg}, {newgen});
%             sgnew.insert({newsg}, {newgen});
%         end
%     end
% end