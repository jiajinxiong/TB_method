classdef Bisection < handle
    properties
        Graph;
        Distance;
        nlayer;
        EdG;Idlayer;
        Lead;
    end

    methods
        function obj = Bisection(G,Lead)
            arguments
                G       graph;
                Lead    (1,:);
            end
            obj.Graph = G;
            obj.Distance = min(distances(G,Lead,'Method','unweighted'));
            obj.nlayer = max(obj.Distance)-1;
            obj.Idlayer = obj.nlayer - obj.Distance+2;
            obj.Lead = Lead;
        end
        function section = partition_graph(obj)
            arguments
                obj     TB_Hamilton.GreenFunction.Bisection;
            end
            LG0 = reshape(string(obj.Graph.Nodes.Name(obj.Distance==max(obj.Distance))),1,[]);
            RG0 = obj.Lead;
            G = rmnode(obj.Graph,[LG0,RG0]);
            [G12,N12,LRG] = obj.bisection(G,obj.nlayer,LG0,RG0);
            PG = G12;     NG = N12;
            while true
                PG_ = {};   N12_ = {};  LRG_ = {};
                for j1 = 1:length(PG)
                    G = PG{j1};         Nlayer = NG{j1};    
                    LG = LRG{j1}{1};    RG = LRG{j1}{2};
                    [G12,N12,LRG12] = obj.bisection(G,Nlayer,LG,RG);
                    PG_(end+1:end+length(G12)) = G12;
                    N12_(end+1:end+length(N12)) = N12;
                    LRG_(end+1:end+length(LRG12)) = LRG12;
                end
                PG = PG_;   NG = N12_; LRG = LRG_;
                if length(PG)==obj.nlayer
                    NPG = cellfun(@(x) x.Nodes.Name,PG,'UniformOutput',false);
                    section = [{LG0},NPG(:)',{RG0}];
                    section = cellfun(@(x) findnode(obj.Graph,x)',section,'UniformOutput',false);
                    break;
                end
            end
        end


        function [BG,N12,LRG] = bisection(obj,G,nlayer,LG,RG)
            arguments
                obj         TB_Hamilton.GreenFunction.Bisection;
                G           graph;
                nlayer      double;
                LG          (1,:) ;
                RG          (1,:) ;
            end
            if nlayer==1
                BG = {G};
                LRG = {[]};
                N12 = {1};
                return;
            end

            node = reshape(string(G.Nodes.Name),1,[]);                  % LG-G-RG
            node_ = [LG,node,RG];
            G0 = subgraph(obj.Graph,node_);
            node0 = G0.Nodes.Name;
            d1 = min(distances(G0,LG),[],1);
            d2 = min(distances(G0,RG),[],1);
            N1 = fix(nlayer/2);     N2 = nlayer - N1;
            N12 = {N1,N2};
            Id1 = d1 <= N1 & d1>0;  Id2 = d2 <= N2 & d2>0;
            node1 = reshape(node0(Id1),1,[]);  node2 = reshape(node0(Id2),1,[]);
            nnodeG1 = G.numnodes;
            Fnode = reshape(setdiff(node,[node1,node2]),1,[]);
            Edge = sort(string(G.Edges.EndNodes),2);

            if length(node1)+length(Fnode) <= N1/nlayer * nnodeG1
                node1 = [node1,Fnode];
                G11 = subgraph(G,node1);   G12 = subgraph(G,node2);
                edge1 = sort(string(G11.Edges.EndNodes),2);
                edge2 = sort(string(G12.Edges.EndNodes),2);
                cut_node = unique(setdiff(Edge,[edge1;edge2],'rows'));
                LR1 = {LG,intersect(cut_node,node2)};
                LR2 = {intersect(cut_node,node1),RG};
                LRG = {LR1,LR2};    BG = {G11,G12};
                return;
            end
            if length(node2)+length(Fnode) <= N2/nlayer * nnodeG1
                node2 = [node2,Fnode];
                G11 = subgraph(G,node1);   G12 = subgraph(G,node2);
                edge1 = sort(string(G11.Edges.EndNodes),2);
                edge2 = sort(string(G12.Edges.EndNodes),2);
                cut_node = unique(setdiff(Edge,[edge1;edge2],'rows'));
                LR1 = {LG,intersect(cut_node,node2)};
                LR2 = {intersect(cut_node,node1),RG};
                LRG = {LR1,LR2};BG = {G11,G12};
                return;
            end
            [G11,G12] = TB_Hamilton.GreenFunction.Bisection.BisectionByKL(G,node1,node2,N1,N2);
            for j1 = 1:10
                [G11_,G12_] = TB_Hamilton.GreenFunction.Bisection.BisectionByKL(G,node1,node2,N1,N2);
                if abs(G11_.numnodes/G12_.numnodes-N1/N2) < abs(G11.numnodes/G12.numnodes-N1/N2)
                    G11 = G11_;  G12 = G12_;
                end
            end
            edge1 = sort(string(G11.Edges.EndNodes),2);
            edge2 = sort(string(G12.Edges.EndNodes),2);
            cut_node = unique(setdiff(Edge,[edge1;edge2],'rows'));
            node1 = G11.Nodes.Name; node2 = G12.Nodes.Name;
            LR1 = {LG,intersect(cut_node,node2)};
            LR2 = {intersect(cut_node,node1),RG};
            LRG = {LR1,LR2};
            BG = {G11,G12};
        end
    end
    methods (Static)

        function [G1,G2] = BisectionByKL(G,node1,node2,N1,N2)
            % Fnode = G - G1 - G2       G = G1--Fnode--G2;
            arguments
                G   graph;
                node1  ;
                node2  ;
                N1  (1,1);
                N2  (1,1);
            end
            Fnode = string(reshape(setdiff(G.Nodes.Name,[node1,node2]),1,[]));
            if isempty(Fnode)
                G1 = subgraph(G,node1);   G2 = subgraph(G,node2);
            else
                numG1 = round(G.numnodes/(N1+N2)*N1);
                % numG2 = G.numnodes - numG1;
                id = randperm(length(Fnode));
                Fnode1 = Fnode(id(1:numG1-length(node1)));
                Fnode2 = Fnode(id(numG1-length(node1)+1:end));
                n = min([length(Fnode1),length(Fnode2)]);
                X = strings(1,n); Y = X;  Gain = zeros(1,n);
                
                Node1 = [node1,Fnode1];             Node2 = [node2,Fnode2];
                G1 = subgraph(G,Node1);             G2 = subgraph(G,Node2);
                for j1 = 1:n
                    g = TB_Hamilton.GreenFunction.Bisection.gain_function(G,G1,G2,Fnode1,Fnode2);
                    [Maxg,id] = max(g,[],'all');
                    Gain(j1) = Maxg;
                    [row,col] = ind2sub(size(g),id);
                    X(j1) = Fnode1(row);  Y(j1) = Fnode2(col);
                    G1 = G1.rmnode(Fnode1(row));
                    G2 = G2.rmnode(Fnode2(col));
                    Fnode1(row) = [];     Fnode2(col) = [];
                end
                sumGain = cumsum(Gain);
                [~,Id] = max(sumGain);
                G1 = subgraph(G,[node1,Fnode1,Y(1:Id),X(Id+1:end)]);  
                G2 = subgraph(G,[node2,Fnode2,X(1:Id),Y(Id+1:end)]);
            end
        end


        function g = gain_function(G,G1,G2,Fnode1,Fnode2)
            arguments
                G       graph;
                G1      graph;
                G2      graph;
                Fnode1  (1,:);
                Fnode2  (1,:);
            end
            D1 = degree(G,Fnode1) - 2*degree(G1,Fnode1);
            D2 = degree(G,Fnode2) - 2*degree(G2,Fnode2);
            C12 = cell2mat(cellfun(@(x) edgecount(G,x,Fnode2)',Fnode1,'UniformOutput',false)');
            g = D1'+D2-C12;
        end
    end
end