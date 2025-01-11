function isGSke = find_W_skeletonic(G)
    arguments
        G digraph;
    end
    id = find(G.Edges.Arrow == 0);
    nid = length(id);
    if nid >=2
        for j1 = 1:nid
            G1 = rmedge(G,id(j1));
            isGSke1 = isscalar(unique(conncomp(G1,'Type','weak')));
            if isGSke1
                for j2 = j1+1:nid
                    G2 = rmedge(G,[id(j1),id(j2)]);
                    isGSke = isscalar(unique(conncomp(G2,'Type','weak')));
                    if ~isGSke
                        return
                    end
                end
            else
                isGSke = isGSke1;
                return
            end
        end
        
    else
        isGSke = true;
    end
end