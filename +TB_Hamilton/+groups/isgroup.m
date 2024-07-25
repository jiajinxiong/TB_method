function Isgroup = isgroup(group)
arguments
    group TB_Hamilton.groups.PointGroupElement;
end
nG = length(group);
G = dictionary(group,1:nG);

Isgroup = arrayfun(@(y) arrayfun(@(x) isKey(G,y*x),group),group,'UniformOutput',false);
Isgroup = all([Isgroup{:}]);
end