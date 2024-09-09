function result = uneq_EE_FG(in_point,out_point,options)
    arguments
        in_point (1,:) string
        out_point (1,:) string = [];
        options.isconnection (1,1) logical = false;
        options.isloop (1,1) logical = true;
        options.isSkeleton (1,1) {mustBeMember(options.isSkeleton,["G","W",""])} = "";
    end
    option = namedargs2cell(options);

    FG = reshape(TB_Hamilton.Feymann_graph.find_all_ee_FG(in_point,out_point,option{:}),[],1);
    FG_add_path = cellfun(@add_path,FG,'UniformOutput',false);
    FG_add_path = cellfun(@(x) TB_Hamilton.Feymann_graph.Feymann_graph(x),FG_add_path,'UniformOutput',false);
    
    Dic = dictionary(FG_add_path{1},1);
    id = 1;
    for j1 = 1:length(FG_add_path)
        if ~isKey(Dic,FG_add_path{j1})
            Dic = insert(Dic,FG_add_path{j1},1);
            id = [id,j1];
        end
    end
    
    
    G1 = Dic.keys;
    num = zeros(1,length(G1));
    for j1 = 1:length(G1)
        num(j1) = nnz(G1(j1)==FG_add_path);
    end
    
    
    
    result = table(FG(id),num','VariableNames',{'uneq_Feymann_Graph','Number'});
    end
    
    
    function G = add_path(G1)
    edges = G1.Edges;
    id = edges.LineStyle == "--";
    G = G1.addedge(edges.EndNodes(id,2),edges.EndNodes(id,1));
    
    end
    