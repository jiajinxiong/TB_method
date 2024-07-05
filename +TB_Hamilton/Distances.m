function d = Distance(G,H0_tag)
    arguments
        G           graph;
        H0_tag      (1,:);
    end
    d = min(distances(G,H0_tag,'Method','unweighted'));
end