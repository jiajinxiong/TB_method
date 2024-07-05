classdef Green_Function < handle
    properties
        Hamilton
        system_graph
        Leads
        GF_pos
        Distance
        PartitionMethod
        MSelfEnergy
    end
    properties (Hidden)
        MItera_Gi_ii;
        MGreenFunction
    end
    methods
        function obj = Green_Function(syst,method)
            arguments
                syst            (1,1)   TB_Hamilton.Builder;
                method   (1,:) char {mustBeMember(method,{'Distance','FM'})} = 'Distance';
            end
            clearAllMemoizedCaches;
            obj.Hamilton = syst.Hamilton;
            obj.system_graph = syst.system_graph;
            obj.Leads = syst.Leads;
            obj.GF_pos = syst.ham_pos;
            obj.PartitionMethod = method;
            obj.sort_layer;
            obj.MItera_Gi_ii = memoize(@obj.Itera_Gi_ii);
            obj.MItera_Gi_ii.CacheSize = max(obj.Distance)+1;
            obj.MGreenFunction = memoize(@obj.GreenFunction);
        end

        function sort_layer(obj)
            if isempty(obj.Leads)
                site_tag = string(obj.system_graph.Nodes.Name);
                site_vec = TB_Hamilton.Tag.get_vec_tag(site_tag);
                ind = (site_vec(1,:) == min(site_vec(1,:)));
                H0_tag = site_tag(ind);
                PG = TB_Hamilton.GreenFunction.Bisection(obj.system_graph,H0_tag);
            else
                PG = TB_Hamilton.GreenFunction.Bisection(obj.system_graph,[obj.Leads.node_H0]);
            end
            if obj.PartitionMethod == "FM"
                D = PG.partition_graph;
                obj.Hamilton = obj.Hamilton([D{:}],[D{:}]);
                obj.system_graph = reordernodes(obj.system_graph,[D{:}]);
                Dis = num2cell(0:length(D)-1);
                obj.Distance = cellfun(@(x,y) repmat(y,1,length(x)),D,Dis,'UniformOutput',false);
                obj.Distance = [obj.Distance{:}];
                obj.GF_pos = obj.GF_pos([D{:}],:);
            else
                D = PG.Idlayer;
                [obj.Distance,ind] = sort(D);
                obj.Hamilton = obj.Hamilton(ind,ind);
                obj.system_graph = reordernodes(obj.system_graph,ind);
                obj.GF_pos = obj.GF_pos(ind,:);
            end

        end

        function [node_gr,gr] = get_lead_gr(obj,lead,E)
            arguments
                obj     TB_Hamilton.GreenFunction.Green_Function;
                lead    struct;
                E       ;
            end
            node0 = [lead.node_H0];
            node1 = [lead.node_V01];
            Inode0 = findnode(obj.system_graph,node0);
            Inode1 = findnode(obj.system_graph,node1);
            H0 = obj.Hamilton(Inode0,Inode0);
            V01 = obj.Hamilton(Inode0,Inode1);
            gr = TB_Hamilton.GreenFunction.SurfaceGreenFunction(H0,V01,E);
            node_gr = node0;
        end

        function self_energy = get_SelfEnergy(obj,lead,E)
            Inode0 = findnode(obj.system_graph,[lead.node_H0]);
            nH = obj.system_graph.numnodes;
            self_energy = spalloc(nH,nH,length(Inode0)^2);
            [~,gr] = obj.get_lead_gr(lead,E);
            node1 = [lead.node_V01];
            Inode1 = findnode(obj.system_graph,node1);
            V01 = obj.Hamilton(Inode0,Inode1);
            self_energy(Inode0,Inode0) = V01'*gr*V01;
        end

        function layer_ind = get_layer(obj,num)
            layer_ind = find(obj.Distance == num);
        end

        function G = GreenFunction(obj,E,options)
            arguments
                obj             TB_Hamilton.GreenFunction.Green_Function;
                E               (1,1) double;
                options.method  (1,1)    = 1;
            end
            nlayer = max(obj.Distance);
            E = E+1i*1e-7;      nH = obj.system_graph.numnodes;

            if options.method == 1
                Fn = obj.Itera_Gi_ii(nlayer,E);
                Inode2 = obj.get_layer(nlayer);
                G = spalloc(nH,nH,length(Inode2)^2);
                G(Inode2,Inode2) = Fn;
            elseif options.method == 2
                Fn = obj.MItera_Gi_ii(nlayer,E,"method",1);
                Inode2 = obj.get_layer(nlayer); Inode0 = Inode2;
                t = tabulate(obj.Distance);
                G = spalloc(nH,nH,sum(t(:,2).^2));
                G(Inode2,Inode2) = Fn;      G0 = Fn;
                for j1 = nlayer-1:-1:0
                    Inode1 = Inode2;  G1 = G0;
                    Fn = obj.MItera_Gi_ii(j1,E,"method",1);
                    Inode2 = obj.get_layer(j1);
                    V01 = obj.Hamilton(Inode2,Inode1);
                    G0 = -Fn*V01*G1;
                    G(Inode2,Inode0) = G0;
                end
            else
                Gn = obj.MItera_Gi_ii(nlayer,E,"method",1);
                Inode2 = obj.get_layer(nlayer);
                G(Inode2,Inode2) = Gn;
                for j1 = nlayer-1:-1:0
                    Inode1 = Inode2;
                    Inode2 = obj.get_layer(j1);
                    G0_nn = obj.MItera_Gi_ii(j1,E,"method",1);
                    V01 = obj.Hamilton(Inode2,Inode1);
                    Gn_1 = G0_nn + G0_nn * V01 * Gn * V01'*G0_nn;
                    G(Inode2,Inode2) = Gn_1;    Gn = Gn_1;
                end
            end
        end

        function T = Transmission(obj,nlead1,nlead2,E)
            arguments
                obj         TB_Hamilton.GreenFunction.Green_Function;
                nlead1      double;
                nlead2      double;
                E           double;
            end
            lead1 = obj.Leads(nlead1);  lead2 = obj.Leads(nlead2);
            Inode1 = findnode(obj.system_graph,[lead1.node_H0]);
            Inode2 = findnode(obj.system_graph,[lead2.node_H0]);
            nH = obj.system_graph.numnodes;
            G = obj.MGreenFunction(E);
            self1 = spalloc(nH,nH,length(Inode1)^2);    self2 = self1;
            self1(Inode1,Inode1) = obj.MSelfEnergy(Inode1,Inode1);
            self2(Inode2,Inode2) = obj.MSelfEnergy(Inode2,Inode2);
            Gamma1 = 1i*(self1-self1');
            Gamma2 = 1i*(self2-self2');
            T = trace(Gamma1 * G * Gamma2 * G');
        end

        function result = Current(obj,E)
            G = obj.GreenFunction(E,'method',2);
            Sigma = obj.MSelfEnergy;
            Gamma = 1i*(Sigma-Sigma');  Sigma_in = Gamma;
            Gn = G * Sigma_in * G';
            [X1,X2] = meshgrid(obj.GF_pos(:,1),obj.GF_pos(:,1));
            [Y1,Y2] = meshgrid(obj.GF_pos(:,2),obj.GF_pos(:,2));
            X = X1-X2;  Y = Y1-Y2;
            J_mn = imag(obj.Hamilton.*Gn);
            result = [diag(J_mn * double(X')) diag(J_mn * double(Y'))];
        end
    end
    methods(Hidden)
        function MSE = Cal_SelfEnergy(obj,E)
            nH = obj.system_graph.numnodes;
            MSE = sparse(nH,nH);
            for j1 = 1:length(obj.Leads)
                MSE = MSE + obj.get_SelfEnergy(obj.Leads(j1),E);
            end
        end


        function F = Itera_Gi_ii(obj,num,E,options)
            arguments
                obj             TB_Hamilton.GreenFunction.Green_Function;
                num             (1,1)   double;
                E               (1,1)   double;
                options.method   = [];
            end

            nlayer = max(obj.Distance);
            if num == 0
                Inode0 = obj.get_layer(0);
                F = (E*speye(length(Inode0))-obj.Hamilton(Inode0,Inode0))\speye(length(Inode0));
            elseif num == 1
                Inode1 = obj.get_layer(1);
                if isempty(options.method)
                    F0 = obj.Itera_Gi_ii(0,E);
                else
                    F0 = obj.MItera_Gi_ii(0,E);
                end
                Inode0 = obj.get_layer(0);
                V01 = obj.Hamilton(Inode0,Inode1);
                F = (E*speye(length(Inode1))-obj.Hamilton(Inode1,Inode1)-V01'*F0*V01)\speye(length(Inode1));
            elseif num == nlayer
                Inode2 = obj.get_layer(num);
                if isempty(options.method)
                    F0 = obj.Itera_Gi_ii(num-1,E);
                else
                    F0 = obj.MItera_Gi_ii(num-1,E);
                end
                Inode1 = obj.get_layer(num-1);
                V01 = obj.Hamilton(Inode1,Inode2);
                if isempty(obj.Leads)
                    self_energy = 0;
                else
                    obj.MSelfEnergy = fetchOutputs(parfeval(backgroundPool,...
                            @obj.Cal_SelfEnergy,1,E));
                    self_energy = obj.MSelfEnergy(Inode2,Inode2);
                end
                d = E*speye(length(Inode2))-obj.Hamilton(Inode2,Inode2)-self_energy;
                F = (d-V01'*F0*V01)\speye(length(Inode2));
            else
                Inode2 = obj.get_layer(num);
                if isempty(options.method)
                    F0 = obj.Itera_Gi_ii(num-1,E);
                else
                    F0 = obj.MItera_Gi_ii(num-1,E);
                end
                Inode1 = obj.get_layer(num-1);
                V01 = obj.Hamilton(Inode1,Inode2);
                F = (E*speye(length(Inode2))-obj.Hamilton(Inode2,Inode2)-V01'*F0*V01)\speye(length(Inode2));
            end
        end
    end
end
