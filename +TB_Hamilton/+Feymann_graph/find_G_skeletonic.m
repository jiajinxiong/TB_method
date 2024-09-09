function isGSke = find_G_skeletonic(G)
    arguments
        G digraph;
    end
    id = find(G.Edges.Arrow ~= 0);
    nid = length(id);
    for j1 = 1:nid
        G1 = rmedge(G,id(j1));
        if isscalar(unique(conncomp(G1,'Type','weak')))
            for j2 = j1:nid
                G2 = rmedge(G,[id(j1),id(j2)]);
                if ~isscalar(unique(conncomp(G2,'Type','weak')))
                    isGSke = 0;
                    return
                end
            end
        else
            isGSke = 0;
            return
        end
    end
    isGSke = 1;
end