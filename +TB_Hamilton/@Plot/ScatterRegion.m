function ScatterRegion(syst,options1,options2)
    arguments
        syst TB_Hamilton.Builder
        % Pcell                                   (1,1)   = 0;
        options1.Pcell                           (1,1) = 0;
        % options1.color                         (1,1) string ='k';
        options1.linewidth                       (1,1) double = 1;
        options2.?matlab.graphics.chart.primitive.Scatter;
    end
    options2 = get_options_site(options2);
    if isempty(syst.TranslationSymmetry)
        alp= 1;
    else
        alp = .4;
    end
    if ~isempty(syst.hopping_pos)
        dim = size(syst.hopping_pos,2)/2;
        pos1 = syst.hopping_pos(:,1:dim)';
        pos2 = syst.hopping_pos(:,dim+1:end)';
        pos = get_plot_pos(pos1,pos2);
        patch(pos(:,1),pos(:,2),'k','EdgeAlpha',alp);
        hold on
        if ~isempty(syst.TranslationSymmetry)
            V = adjacency(syst.system_graph);
            [x,y] = gplot(V,syst.ham_pos(syst.cell_ind,:));
            plot(x,y,'k','LineWidth',1);
        end
    end
    sites = syst.Sites;
    site_pos = reshape([sites.pos],2,[])';
    site_vec = TB_Hamilton.Tag.get_vec_tag(arrayfun(@(x) x.tag(1),sites))';
    hold on;

    for j2 = 1:length(syst.site_family)
        ind = [sites.site_family]==syst.site_family(j2);
        site_pos1 = site_pos(ind,:);
        if length(options2)<=1
            propertyCell = namedargs2cell(options2);
        else
            propertyCell = namedargs2cell(options2(j2));
        end
        propertyCell(end+1:end+4) = {'MarkerFaceAlpha',sqrt(alp/100)*10,'MarkerEdgeAlpha',sqrt(alp/100)*10};
        Plot_Sites(site_pos1,site_vec(ind,:),syst.site_family(j2),propertyCell{:});
    end
    if ~isempty(syst.TranslationSymmetry)
        pos1 = syst.ham_pos(syst.cell_ind,1);
        pos2 = syst.ham_pos(syst.cell_ind,2);
        scatter(pos1,pos2,'k','filled');
    end


    syst.hopping_pos=[]; syst.Leads_hopping = [];
    if options1.Pcell
        plot_cell(syst);
    end
    if ~isempty(syst.Leads)
        Plot_Leads(syst);
    end
    axis equal;
    hold off
end

function Plot_Sites(site_pos,site_vec,site_family,options)
    arguments
        site_pos    (:,2)   double;
        site_vec    (:,2)   double;
        site_family;
        options.?matlab.graphics.chart.primitive.Scatter
    end
    propertyCell = namedargs2cell(options);
    s = scatter(site_pos(:,1),site_pos(:,2),'filled',propertyCell{:});
    dtt = s.DataTipTemplate;
    site_family_ = dataTipTextRow('Site Family',repmat(site_family,1,size(site_pos,1)));
    site_vec_1 = dataTipTextRow('a_1',site_vec(:,1)');
    site_vec_2 = dataTipTextRow('a_2',site_vec(:,2)');
    dtt.DataTipRows(1) = site_family_;
    dtt.DataTipRows(2) = site_vec_1;
    dtt.DataTipRows(end+1) = site_vec_2;
end


function Plot_Leads(syst)
    arguments
        syst        TB_Hamilton.Builder;
    end
    dim = size(syst.hopping_pos,2)/2;sites = syst.Sites;
    site_pos = reshape([sites.pos],2,[])';
    for j1 = 1:length(syst.Leads_hopping)
        lead_pos1 = syst.Leads_hopping(j1).pos(:,1:dim)';
        lead_pos2 = syst.Leads_hopping(j1).pos(:,dim+1:end)';
        lead_pos = get_plot_pos(lead_pos1,lead_pos2);
        plot(lead_pos(:,1),lead_pos(:,2),'r','LineWidth',1);
        Lead_site_pos = unique(lead_pos,'rows');
        Lead_site_pos1 = setdiff(Lead_site_pos,site_pos,'rows');
        scatter(Lead_site_pos1(:,1),Lead_site_pos1(:,2),'filled','MarkerFaceColor','r','MarkerFaceAlpha',1);
        x0 = syst.Leads(j1).TranslationSymmetry * syst.vec;
        patch(lead_pos(:,1)+x0(1),lead_pos(:,2)+x0(2),'red','EdgeColor','red','EdgeAlpha',.4,'LineWidth',1);
        scatter(Lead_site_pos(:,1)+x0(1),Lead_site_pos(:,2)+x0(2),'filled','MarkerFaceColor','r',...
            'MarkerFaceAlpha',.4);
    end
end


function plot_cell(syst)
    arguments
        syst    TB_Hamilton.Builder;
    end
    site = mean(double(syst.basises));
    vec = syst.vec;
    vec_tag= TB_Hamilton.Tag.get_vec_tag([syst.Sites.tag]);
    nx = min(vec_tag(1,:)):max(vec_tag(1,:));    
    ny = min(vec_tag(2,:)):max(vec_tag(2,:));
    [XX,YY] = meshgrid(nx,ny);
    Site = double(site + [XX(:),YY(:)] * vec);
    [Vx,Vy] = voronoi(Site(:,1),Site(:,2));
    Vx(abs(Vx)>max(abs(Site(:,1)))) = NaN;
    Vy(abs(Vy)>max(abs(Site(:,2)))) = NaN;
    plot(Vx,Vy,'--','color',[.5,.5,.5]);
end


function result = get_options_site(option)
    if ~isempty(fieldnames(option))
        [~,ia] = max(structfun(@(x) length(x),option));
        fields = fieldnames(option);
        A = option.(fields{ia});
        result = struct(fields{ia},num2cell(A));
        ind_ = setdiff(1:length(fields),ia);
        for j = ind_
            A = option.(fields{j});
            if length(A)==1 || isempty(A)
                result = arrayfun(@(x) setfield(x,fields{j},A),result);
            else
                result = arrayfun(@(x,y) setfield(x,fields{j},y),result,A);
            end
        end
    else
        result = struct();
    end
end

function pos_new = get_plot_pos(pos1,pos2)
    arguments
        pos1    (2,:)   double;
        pos2    (2,:)   double;
    end
    posx = reshape([pos1(1,:);pos2(1,:)],[],1);
    posy = reshape([pos1(2,:);pos2(2,:)],[],1);
    pos_new = zeros(3*length(pos1),2);
    ind_NaN = 3:3:length(pos_new);
    pos_new(ind_NaN,:) = nan(length(ind_NaN),2);
    pos_new(~isnan(pos_new(:,1)),:) = [posx posy];
end