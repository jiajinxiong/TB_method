function FG = find_all_ep_FG(in_point,out_point,options)
arguments
    in_point (1,:) string
    out_point (1,:) string = [];
    options.isconnection (1,1) logical = false;
    options.isloop (1,1) logical = true;
    options.isGSke (1,1) logical = false;
end

isconnection = options.isconnection;
isloop = options.isloop;
isGSke = options.isGSke;
in_ele = in_point(contains(in_point,'_e'));
in_p = setdiff(in_point,in_ele);
out_ele = out_point(contains(out_point,'e'));
out_p = setdiff(out_point,out_ele);
all_ele = [in_ele,out_ele];     % all electron
all_p = [in_p,out_p];           % all phonon

ce = [in_ele(contains(in_ele,'c')),out_ele(contains(out_ele,'i'))];  % create the electron
ae = setdiff(all_ele,ce);       % annihilation electron
cp = [in_p(contains(in_p,'c')),out_p(contains(out_p,'i'))];  % create the phonon
ap = setdiff(all_p,cp);          % annihilation phonon

Paes = perms(ae);           % permutation of annihilation electron
Paps = perms(ap);           % permutation of annihilation phonon

rPaes = regexprep(Paes,'_.*','');   % remove the suffix of the permutation of annihilation electron
rPaps = regexprep(Paps,'_.*','');   % remove the suffix of the permutation of annihilation phonon
rce = regexprep(ce,'_.*','');
rcp = regexprep(cp,'_.*','');


% if ~isloop
%     id = cellfun(@(x) ~any([rce,rcp] == string(x.Edges.EndNodes(:,2))),FG);
%     FG = FG(id);
% end


nrPaes = size(rPaes,1);
nrPaps = size(rPaps,1);
FG = cell(nrPaps,nrPaes);

Edge_Style_Template = repmat("-", 1, numel(rce) + numel(rcp));
Arrow_Template = repmat(12, 1, numel(rce) + numel(rcp));


parfor j1 = 1:nrPaps
    localFG = cell(1, nrPaes);
    for j2 = 1:nrPaes
        G1 = digraph([rce,rcp],[rPaes(j2,:),rPaps(j1,:)]);
        [idxOut,m] = findedge(G1,rcp,rPaps(j1,:));
        [~,id] = unique(m);
        pp = idxOut(id);

        Edge_Style = Edge_Style_Template;
        Edge_Style(pp) = "--";
        Arrow = Arrow_Template; Arrow(pp) = 0;
        
        G1.Edges.LineStyle = Edge_Style';
        G1.Edges.Arrow = Arrow';
        localFG{j2} = G1;
    end
    FG(j1, :) = localFG;
end
if isconnection
    id = cellfun(@(x) isscalar(unique(conncomp(x,'Type','weak'))),FG);
    FG = FG(id);
end

% This step can be futher optimized.
if ~isloop
    id = cellfun(@(x) ~any(string(x.Edges.EndNodes(:,1)) == string(x.Edges.EndNodes(:,2))),FG);
    FG = FG(id);
end

if isGSke
    id = cellfun(@(x) TB_Hamilton.Feymann_graph.find_G_skeletonic(x),FG);
    FG = FG(find(id));
end


end