function is_Irr = find_irreducible(G,tag,value)
    arguments
        G digraph;
        tag string;
        value ;
    end
    id = find(G.Edges.(tag) == value);
    nid = length(id);
    for j1 = 1:nid
        G1 = rmedge(G,id(j1));
        if ~isscalar(unique(conncomp(G1,'Type','weak')))
            is_Irr = 0;
            return
        end
    end
    is_Irr = 1;
end