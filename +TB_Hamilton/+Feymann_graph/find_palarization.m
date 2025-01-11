function P = find_palarization(G)
    arguments
        G digraph;
    end
    P = {};
    id = find(G.Edges.Arrow == 0);
    nid = length(id);
    Names = "ep"+(1:10);
    for j1 = 1:nid
        % G1 = rmedge(G,id(j1));
        for j2 = j1+1:nid
            
            G2 = rmedge(G,[id(j1),id(j2)]);
            G2_conncomp = conncomp(G2,'Type','weak');
            if ~isscalar(unique(G2_conncomp))
                G2_1 = subgraph(G2,G2_conncomp == 1);
                G2_2 = subgraph(G2,G2_conncomp == 2);
                degree_G2_1 = indegree(G2_1) + outdegree(G2_1);
                degree_G2_2 = indegree(G2_2) + outdegree(G2_2);
                id1 = degree_G2_1 == 2; id2 = degree_G2_2 == 2;
                G2_1.Nodes.Name(id1) = {'ei','ef'};
                G2_2.Nodes.Name(id2) = {'ei','ef'};
                is_Irr_G2_1 = TB_Hamilton.Feymann_graph.find_irreducible(G2_1,'Arrow',0);
                is_Irr_G2_2 = TB_Hamilton.Feymann_graph.find_irreducible(G2_2,'Arrow',0);
                n1 = G2_1.numnodes; n2 = G2_2.numnodes;
                G2_1.Nodes.Name(~id1) = cellstr(Names(1:n1-2));
                G2_2.Nodes.Name(~id2) = cellstr(Names(1:n2-2));

                if is_Irr_G2_1
                    P{end+1}=G2_1;
                end
                if is_Irr_G2_2
                    P{end+1}=G2_2;
                end
            end
        end
    end
end