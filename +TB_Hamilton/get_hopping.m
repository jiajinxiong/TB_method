function [tag,value,pos,sites] = get_hopping(hopping_site, hopping_value)
    arguments
        hopping_site            struct;
        hopping_value           cell;
    end
    size_hopping = size(hopping_site);
    if isa(hopping_value{1},'cell')
        hopping_value = hopping_value{1};
    elseif length(hopping_value) ~= size_hopping(1)
        hopping_value(1:size_hopping(1)) = hopping_value;
    end
    if size_hopping(2) == 1
        [~,ia] = unique([hopping_site.site_family]);
        A_ = {hopping_site(ia).tag};
        orbits = cellfun(@(x) length(x),A_);
        pos = [];
        sites = hopping_site;
        if max(orbits) == 1
            value = cell2mat(hopping_value);
            tag = [[hopping_site.tag] ;[hopping_site.tag]];
        else
            value = cell2mat(cellfun(@(x) reshape(x,1,[]),...
                hopping_value,'UniformOutput',false));
            A = {hopping_site.tag};
            if min(orbits)==max(orbits)
                A1 = cellfun(@(x) reshape(repmat(x',1,orbits(1))',1,[]),A,'UniformOutput',false);
            else
                A_orbits = cellfun(@(x) length(x),A,'UniformOutput',false);
                A1 = cellfun(@(x,y) reshape(repmat(x',1,y)',1,[]),A,A_orbits,'UniformOutput',false);
            end
            A1 = [A1{:}];   ind = value~=0;
            tag = [A1(ind);A1(ind)];
            value = value(ind);
        end
    elseif size_hopping(2) == 2
        hoppig_site_1 = hopping_site(:,1);
        hoppig_site_2 = hopping_site(:,2);
        sites = [hoppig_site_1;hoppig_site_2];
        [~,ia] = unique([sites.site_family]);
        pos1 = cell2mat({hoppig_site_1.pos}');
        pos2 = cell2mat({hoppig_site_2.pos}');
        ind = min(pos1==pos2,[],2);
        pos = [pos1(~ind,:) pos2(~ind,:)];
        [~,ind1] = unique([pos1;pos2],'rows');
        sites = sites(ind1);
        A = {hoppig_site_1.tag};
        B = {hoppig_site_2.tag};
        A_orbits = cellfun(@(x) length(x),A,'UniformOutput',false);
        B_orbits = cellfun(@(x) length(x),B,'UniformOutput',false);
        orbits = [A_orbits,B_orbits];
        orbits = cell2mat(orbits(ia));
        if max(orbits) == 1
            value = cell2mat(hopping_value);ind = value~=0;
            value = value(ind);
            tag = [[hoppig_site_1(ind).tag] ;[hoppig_site_2(ind).tag]];
            
        else
            value = cell2mat(cellfun(@(x) reshape(x,1,[]),...
                hopping_value,'UniformOutput',false))';
            A1 = cellfun(@(x,y) reshape(repmat(x',1,y)',1,[]),A,B_orbits,'UniformOutput',false);
            B1 = cellfun(@(x,y) repmat(x,1,y),B,A_orbits,'UniformOutput',false);
            B1 = [B1{:}];   A1 = [A1{:}];   value = value(:);
            ind = value~=0;
            value = value(ind);
            tag = [A1(ind);B1(ind)];
        end
    else
        error("Hopping Site's size must be (:,2) or (:,1)")
    end
end