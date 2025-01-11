function FG = find_all_ee_FG(in_point,out_point,options)
    arguments
        in_point (1,:) string
        out_point (1,:) string = [];
        options.isconnection (1,1) logical = false;
        options.isloop (1,1) logical = true;
        options.isSkeleton (1,1) {mustBeMember(options.isSkeleton,["G","W","","GW"])} = "";
    end
    isconnection = options.isconnection;
    isloop = options.isloop;
    isGSke = options.isSkeleton;

    in_ele = in_point(contains(in_point,'_e'));
    in_p = setdiff(in_point,in_ele);
    out_ele = out_point(contains(out_point,'e'));
    out_p = setdiff(out_point,out_ele);
    all_ele = [in_ele,out_ele];     % all electron
    all_p = [in_p,out_p];           % all phonon
    
    ce = [in_ele(contains(in_ele,'c')),out_ele(contains(out_ele,'i'))];  % create the electron
    ae = setdiff(all_ele,ce);       % annihilation electron
    cp = [in_p(contains(in_p,'c')),out_p(contains(out_p,'i'))];  % create the phonon
    ap = setdiff(all_p,cp);          % annihilation photon
    
    % Generate permutations of annihilation electrons
    Paes = perms(ae);           % permutation of annihilation electron
    % Paps = perms(ap);           % permutation of annihilation photon
    
    % Remove suffixes from particle identifiers
    rPaes = regexprep(Paes,'_.*','');   % remove the suffix of the permutation of annihilation electron
    rap = regexprep(ap,'_.*','');   % remove the suffix of the permutation of annihilation photon
    rce = regexprep(ce,'_.*','');
    rcp = regexprep(cp,'_.*','');
    
    % Initialize variables for graph construction
    nrPaes = size(rPaes,1);
    FG = cell(1, nrPaes);
    % if isGSke == "G"
    %     Edge_Style_Template = repmat("G-", 1, numel(rce) + numel(rcp));
    % else
    %     Edge_Style_Template = repmat("-", 1, numel(rce) + numel(rcp));
    % end
    Arrow_Template = repmat(20, 1, numel(rce) + numel(rcp));
    Edge_Style_Template = repmat("-", 1, numel(rce) + numel(rcp));

    
    % Parallel loop to construct directed graphs
    parfor j1 = 1:nrPaes
        G1 = digraph([rce,rcp],[rPaes(j1,:),rap]);

        nodes = G1.Nodes.Name;
        in_out = strings(size(nodes));
        id = [find(contains(nodes,'i'));find(contains(nodes,'f'))];
        in_out(id) = nodes(id);
        in_out(find(cellfun('isempty', cellstr(in_out)))) = "in";
        G1.Nodes.in_out = in_out;
        
        [idxOut,m] = findedge(G1,rcp,rap);
        [~,id] = unique(m);
        pp = idxOut(id);

        Edge_Style = Edge_Style_Template;
        Edge_Style(pp) = "--";
        Arrow = Arrow_Template; Arrow(pp) = 0;
        
        G1.Edges.LineStyle = Edge_Style';
        G1.Edges.Arrow = Arrow';
        FG{j1} = G1;
    end

    % Filter graphs based on connection option
    if isconnection
        id = cellfun(@(x) isscalar(unique(conncomp(x,'Type','weak'))),FG);
        FG = FG(id);
    end
    
    % Filter graphs based on loop option
    % This step can be futher optimized.
    if ~isloop
        id = cellfun(@(x) ~any(string(x.Edges.EndNodes(:,1)) == string(x.Edges.EndNodes(:,2))),FG);
        FG = FG(id);
    end
    
    if isGSke == "G"
        id = cellfun(@(x) TB_Hamilton.Feymann_graph.find_G_skeletonic(x),FG);
        FG = FG((id));
    elseif isGSke == "W"
        id = cellfun(@(x) TB_Hamilton.Feymann_graph.find_W_skeletonic(x),FG);
        FG = FG(id);
        % FG = cellfun(@(x) TB_Hamilton.Feymann_graph.find_palarization(x),FG,'UniformOutput',false);
        % FG = cat(2, FG{:});
    elseif isGSke == "GW"
        id = cellfun(@(x) TB_Hamilton.Feymann_graph.find_G_skeletonic(x),FG);
        FG = FG(id);
        id = cellfun(@(x) TB_Hamilton.Feymann_graph.find_W_skeletonic(x),FG);
        FG = FG(id);
        
    end
    
end