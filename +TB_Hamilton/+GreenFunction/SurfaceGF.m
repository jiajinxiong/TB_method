classdef SurfaceGF < handle
    %SURFACEGF 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        H0_ind;
        V1_ind;
        Ham;
    end
    
    methods
        function obj = SurfaceGF(syst,Translation)
            arguments
                syst            TB_Hamilton.Builder;
                Translation     (1,:)   double;
            end
            nodes_tag = string(syst.system_graph.Nodes.Name);
            nodes_tag_Trans = TB_Hamilton.Tag.get_tag_TranslationSymmetry(nodes_tag,Translation);
            H0_tag = setdiff(nodes_tag,nodes_tag_Trans);
            V1_tag = TB_Hamilton.Tag.get_tag_TranslationSymmetry(H0_tag,Translation);
            obj.H0_ind = findnode(syst.system_graph,H0_tag);
            obj.V1_ind = findnode(syst.system_graph,V1_tag);
            obj.get_ham(syst);
        end
        
        function get_ham(obj,syst)
            if isa(syst.Hamilton,"double")
                obj.Ham = syst.Hamilton;
            elseif isa(syst.Hamilton,"function_handle")
                obj.Ham = @(k) syst.Hamilton(k);
            end
        end
        function gr = SGF(obj,E,k)
            arguments
                obj     TB_Hamilton.GreenFunction.SurfaceGF;
                E       (1,1)   double;
                k       (1,1)   double = [];
            end
            if isempty(k)
                H0 = obj.Ham(obj.H0_ind,obj.H0_ind);  V1 = obj.Ham(obj.H0_ind,obj.V1_ind);
            elseif length(k)==1
                Ham_k = obj.Ham(k);
                H0 = Ham_k(obj.H0_ind,obj.H0_ind);  V1 = Ham_k(obj.H0_ind,obj.V1_ind);
            end
            gr = TB_Hamilton.GreenFunction.SurfaceGreenFunction(H0,V1,E);
        end
    end
end

