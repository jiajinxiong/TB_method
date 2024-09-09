classdef Feymann_graph

properties
    G;
end

methods
    function obj = Feymann_graph(G)  
        arguments
            G (1,1) digraph;
        end
        obj.G = G;
    end
    function iseq = eq(obj,FG)
        arguments
            obj     TB_Hamilton.Feymann_graph.Feymann_graph;
            FG      ;
        end
        if iscell(FG)
            iseq = cellfun(@(x) isisomorphic(obj.G,x.G),FG);
        else
            iseq = isisomorphic(obj.G,FG.G);
        end
    
    end
    function result = keyMatch(obj,FG)
        arguments
            obj     ;
            FG      ;
        end
        result = obj==FG;
    end
    function result = keyHash(obj)
        if iscell(obj)
            result = keyHash(sort(obj{1}.G.Nodes.Name));
        else
            result = keyHash(sort(obj.G.Nodes.Name));
        end
        
    end
end
end