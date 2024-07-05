classdef Builder < handle
    properties
        site_family
        system_graph
        orbits
        Hamilton
        ham_pos
        Leads
        basises
        vec
        TranslationSymmetry
        cell_ind
    end
    properties (Hidden)
        hopping_tag
        hopping_pos
        Leads_hopping
        hopping_value
        Sites
    end
    methods
        function obj = Builder(Lattice,TranslationSymmetry)
            arguments
                Lattice;
                TranslationSymmetry  = [];
            end
            obj.vec = Lattice.prim_vecs;
            obj.basises = Lattice.basises;
            obj.hopping_value = [];
            obj.hopping_tag = [];
            obj.site_family = Lattice.names;
            obj.system_graph;   obj.Sites;
            obj.orbits = Lattice.orbits; 
            obj.Hamilton;   obj.ham_pos;
            obj.TranslationSymmetry = TranslationSymmetry;
        end

        function obj = subsasgn(obj, S, B)
            if strcmp(S(1).type, '()') && isa(obj,'TB_Hamilton.Builder')
                obj.set_hopping(S.subs{:},B)
            end
        end

        function set_hopping(obj,varargin)
            import TB_Hamilton.get_hopping
            if isa(varargin{end},'function_handle')
                fun = varargin{end};
                if iscell(varargin{1})
                    varargin(end) = {fun(varargin{1}{:})};
                else
                    varargin(end) = {fun(varargin{1:end-1})};
                end
            end

            if iscell(varargin{1})
                    [tag,value,pos,sites] = get_hopping([varargin{1}{1},varargin{1}{2}],varargin(end));
            else
                if nargin==3 && min(size(varargin{1}))==1
                    [tag,value,pos,sites] = get_hopping(varargin{1}, varargin(end));
                elseif min(size(varargin{1}))==2
                    [tag,value,pos,sites] = get_hopping(varargin{1}, varargin(end));
                else
                    [tag,value,pos,sites] = get_hopping([varargin{1} varargin{2}], varargin(end));
                end
            end
            obj.hopping_tag = [obj.hopping_tag tag];
            obj.hopping_value = [obj.hopping_value,reshape(value,1,[])];
            
            size_pos = size(pos);
            dim = size_pos(2)/2;
            if isempty(obj.hopping_pos)
                obj.Sites = sites;
            else
                if ~isempty(pos)
                    pos1 = [obj.hopping_pos(:,1:dim);obj.hopping_pos(:,dim+1:end)];
                    pos2 = cell2mat({sites.pos}');
                    [~,ind] =setdiff(pos2,pos1,'rows');
                    obj.Sites = [obj.Sites;sites(ind)];
                end

            end
            obj.hopping_pos = [obj.hopping_pos;pos];
        end

        function attach_lead(obj,Lead)
            arguments
                obj
                Lead    TB_Hamilton.Builder
            end
            if isempty(obj.Leads)
                obj.Leads_hopping.pos = Lead.hopping_pos;
                obj.Leads.TranslationSymmetry = Lead.TranslationSymmetry;
                obj.Leads.node_H0 = setdiff([Lead.Sites.tag],[obj.Sites.tag]);
                obj.Leads.sites = Lead.Sites;
                obj.Leads.node_V01 = TB_Hamilton.Tag.get_tag_TranslationSymmetry(obj.Leads.node_H0,-Lead.TranslationSymmetry)';
            else
                nLeads = length(obj.Leads);
                obj.Leads_hopping(nLeads+1).pos = Lead.hopping_pos;
                obj.Leads(nLeads+1).TranslationSymmetry = Lead.TranslationSymmetry;
                obj.Leads(nLeads+1).node_H0 = setdiff([Lead.Sites.tag],[obj.Sites.tag]);
                obj.Leads(nLeads+1).sites = Lead.Sites;
                obj.Leads(nLeads+1).node_V01 = TB_Hamilton.Tag.get_tag_TranslationSymmetry(obj.Leads(end).node_H0,-Lead.TranslationSymmetry)';
            end
            obj.hopping_tag = [obj.hopping_tag Lead.hopping_tag];
            obj.hopping_value = [obj.hopping_value,Lead.hopping_value];
        end

        function finalized(obj)
            if ~isempty(obj.Leads)
                leads_sites = {obj.Leads.sites}';
                leads_sites = cell2mat(leads_sites);
                lead_sites_pos = cell2mat({leads_sites.pos}');
                [~,ind] =setdiff(lead_sites_pos,cell2mat({obj.Sites.pos}'),'rows');
                all_sites = [obj.Sites;leads_sites(ind)];
            else
                all_sites = obj.Sites;
            end
            tag = {all_sites.tag};    site_pos = {all_sites.pos}';
            orbit = cellfun(@(x) length(x),tag,'UniformOutput',false);
            pos = cell2mat(cellfun(@(x,y) repmat(x,y,1),site_pos,orbit','UniformOutput',false));
            tag = [tag{:}];
            obj.hopping_value = double(obj.hopping_value);
            if max(abs(imag(obj.hopping_value)))==0
                SG = graph(obj.hopping_tag(1,:),obj.hopping_tag(2,:),obj.hopping_value);
                obj.Hamilton = adjacency(SG,'weighted');
                obj.system_graph = graph(obj.hopping_tag(1,:),obj.hopping_tag(2,:),'omitselfloops');
                ind = findnode(obj.system_graph,tag);
                [~,ind_] = sort(ind);
                obj.ham_pos = pos(ind_,:);
            else
                SG = graph(obj.hopping_tag(1,:),obj.hopping_tag(2,:),real(obj.hopping_value));
                G1 = digraph(obj.hopping_tag(1,:),obj.hopping_tag(2,:),imag(obj.hopping_value));
                H0 = adjacency(SG,'weighted');
                H1 = 1i*adjacency(G1,'weighted');
                obj.Hamilton = H0 + H1+H1';
                obj.system_graph = graph(obj.hopping_tag(1,:),obj.hopping_tag(2,:),'omitselfloops');
                ind = findnode(obj.system_graph,tag);
                [~,ind_] = sort(ind);
                obj.ham_pos = pos(ind_,:);
            end
            Ham_ = obj.Hamilton;

            if size(obj.TranslationSymmetry,1) == 1
                translation_node = TB_Hamilton.Tag.get_tag_TranslationSymmetry([obj.Sites.tag],obj.TranslationSymmetry);
                H0_tag = setdiff([obj.Sites.tag],translation_node);  % the most small cell
                H0_hopping_vec_relative = cellfun(@(x) get_vec_tag(neighbors(obj.system_graph,x)')-get_vec_tag(x),...
                    H0_tag,'UniformOutput',false);
                H0_hopping_vec_relative = unique([H0_hopping_vec_relative{:}]','rows')';
                hopping_cell = unique(obj.TranslationSymmetry * H0_hopping_vec_relative)...
                    /norm(obj.TranslationSymmetry)^2;
                hopping_cell = hopping_cell(hopping_cell>0);
                hopping_cell = unique(ceil(abs(hopping_cell)));
                ind0 = findnode(obj.system_graph,H0_tag);
                obj.cell_ind = ind0;
                ind1 = zeros(length(H0_tag),length(hopping_cell));
                for j1 = 1:length(hopping_cell)
                    H1_tag = TB_Hamilton.Tag.get_tag_TranslationSymmetry(H0_tag,...
                        obj.TranslationSymmetry * hopping_cell(j1));
                    ind1(:,j1) = findnode(obj.system_graph,H1_tag);
                end
                obj.system_graph = subgraph(obj.system_graph,H0_tag);
                obj.Hamilton = @(x) obj.Ham_Ribbon(Ham_,hopping_cell,ind1,kx=x);

            elseif size(obj.TranslationSymmetry,1) == 2
                site_tag = string(table2array(obj.system_graph.Nodes));
                site_vec_tag = get_vec_tag(site_tag);
                mid_site = fix(mean(site_vec_tag,2));
                H0_tag = site_tag(ismember(site_vec_tag',mid_site','rows'));
                H0_hopping_vec_relative = cellfun(@(x) neighbors(obj.system_graph,x)',H0_tag,'UniformOutput',false);
                H0_hopping_vec_relative = setdiff(string([H0_hopping_vec_relative{:}]),H0_tag);
                coordinates = get_vec_tag(H0_hopping_vec_relative);
                hopping_cell = unique((obj.TranslationSymmetry * coordinates)',"rows")-mid_site';
                obj.ham_pos = double(obj.ham_pos);
                ind0 = findnode(obj.system_graph,H0_tag);
                obj.cell_ind = ind0;
                ind1 = zeros(length(H0_tag),length(hopping_cell));
                for j1 = 1:length(hopping_cell)
                    H1_tag = TB_Hamilton.Tag.get_tag_TranslationSymmetry(H0_tag,hopping_cell(j1,:));
                    ind1(:,j1) = findnode(obj.system_graph,H1_tag);
                end
                obj.system_graph = subgraph(obj.system_graph,H0_tag);
                obj.Hamilton = @(kx,ky) obj.Hamilton_bulk(Ham_,hopping_cell,ind1,kx,ky);
            end
            obj.hopping_tag=[];obj.hopping_value=[];    %obj.hopping_pos=[];

            function coordinates = get_vec_tag(tag)
                matches = regexp(tag, '\((.*?)\)', 'tokens'); % 使用正则表达式匹配括号中的内容
                coordinates_str = cellfun(@(x) x{1},matches,'UniformOutput',false); % 获取匹配到的内容
                coordinates = cellfun(@(x) str2num(strrep(x, ',',' '))',string(coordinates_str),'UniformOutput',false);
                coordinates = [coordinates{:}];
            end
        end
        H = Ham_Ribbon(obj,Ham_,hopping_cell,ind1,Mk);

        function H = Hamilton_bulk(obj,Ham_,hopping_cell,ind1,kx,ky)
            ind0 = obj.cell_ind;
            x0 = obj.ham_pos(ind0,1)-obj.ham_pos(ind0,1)';
            y0 = obj.ham_pos(ind0,2)-obj.ham_pos(ind0,2)';
            % pos0 = obj.ham_pos(ind0,:) * [kx;ky];
            Pos = x0*kx+y0*ky;
            % [Pos0,Pos1] = meshgrid(pos0,pos0);
            % Pos = Pos0-Pos1;
            H = Ham_(ind0,ind0).* exp(1i * Pos);
            for j1 = 1:size(hopping_cell,1)
                H1 = Ham_(ind0,ind1(:,j1));
                x0 = obj.ham_pos(ind0,1)-obj.ham_pos(ind1(:,j1),1)';
                y0 = obj.ham_pos(ind0,2)-obj.ham_pos(ind1(:,j1),2)';
                Pos = x0 *kx+y0*ky;
                % pos1 = obj.ham_pos(ind1(:,j1),:) * [kx;ky];
                % [Pos0,Pos1] = meshgrid(pos0,pos1);
                % Pos = Pos0-Pos1;
                H = H + H1 .* exp(1i * Pos);
            end
            H = (H+H')/2;
        end
    end
    methods(Static)
        function hopping = HoppingKind(lat1,lat2,kind)
            arguments
                lat1    TB_Hamilton.Monatomic;
                lat2    TB_Hamilton.Monatomic;
                kind    (1,:)   double;
            end
            lat1_tag = lat1.LatInShape;
            lat2_tag = lat1_tag+kind;
            [bx,by] = boundary(lat2.shape);
            Point2x = lat2_tag(:,1) * lat2.prim_vecs(1,1) + ...
                lat2_tag(:,2) * lat2.prim_vecs(2,1) + lat2.basises(1);
            Point2y = lat2_tag(:,1) * lat2.prim_vecs(1,2) + ...
                lat2_tag(:,2) * lat2.prim_vecs(2,2) + lat2.basises(2);
            in = inpolygon(Point2x,Point2y,bx,by);
            hopping{1} = lat1(lat1_tag(in(:),1),lat1_tag(in(:),2));
            hopping{2} = lat2(lat2_tag(in(:),1),lat2_tag(in(:),2));
        end
    end
end