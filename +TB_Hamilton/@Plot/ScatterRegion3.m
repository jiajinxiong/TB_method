function ScatterRegion3(syst,options1,options2)
arguments
    syst        TB_Hamilton.Builder
    options1.linewidth                       (1,1) double = 1;
    options2.?matlab.graphics.chart.primitive.Scatter;
end
if ~isempty(syst.hopping_pos)
    dim = size(syst.hopping_pos,2)/2;
    pos1 = syst.hopping_pos(:,1:dim)';
    pos2 = syst.hopping_pos(:,dim+1:end)';
    pos = get_plot_pos(pos1,pos2);
    plot3(pos(:,1),pos(:,2),pos(:,3),'k');
    hold on
    % if ~isempty(syst.TranslationSymmetry)
    %     V = adjacency(syst.system_graph);
    %     [x,y] = gplot(V,syst.ham_pos(syst.cell_ind,:));
    %     plot3(x,y,'k','LineWidth',1);
    % end
end
sites = syst.Sites;
site_pos = reshape([sites.pos],3,[])';
site_vec = TB_Hamilton.Tag.get_vec_tag(arrayfun(@(x) x.tag(1),sites))';
for j2 = 1:length(syst.site_family)
    ind = [sites.site_family]==syst.site_family(j2);
    site_pos1 = site_pos(ind,:);
    if length(options2)<=1
        propertyCell = namedargs2cell(options2);
    else
        propertyCell = namedargs2cell(options2(j2));
    end
    % propertyCell(end+1:end+4) = {'MarkerFaceAlpha',sqrt(alp/100)*10,'MarkerEdgeAlpha',sqrt(alp/100)*10};
    Plot_Sites(site_pos1,site_vec(ind,:),syst.site_family(j2),propertyCell{:});
end

end

function Plot_Sites(site_pos,site_vec,site_family,options)
arguments
    site_pos    (:,3)   double;
    site_vec       double;
    site_family;
    options.?matlab.graphics.chart.primitive.Scatter
end
propertyCell = namedargs2cell(options);
s = scatter3(site_pos(:,1),site_pos(:,2),site_pos(:,3),'filled',propertyCell{:});
dtt = s.DataTipTemplate;
site_family_ = dataTipTextRow('Site Family',repmat(site_family,1,size(site_pos,1)));
site_vec_1 = dataTipTextRow('a_1',site_vec(:,1)');
site_vec_2 = dataTipTextRow('a_2',site_vec(:,2)');
dtt.DataTipRows(1) = site_family_;
dtt.DataTipRows(2) = site_vec_1;
dtt.DataTipRows(end) = site_vec_2;
if size(site_vec,2) == 3
    site_vec_3 = dataTipTextRow('a_3',site_vec(:,3)');
    dtt.DataTipRows(end+1) = site_vec_3;
end


end

function pos_new = get_plot_pos(pos1,pos2)
arguments
    pos1    (3,:)   double;
    pos2    (3,:)   double;
end
posx = reshape([pos1(1,:);pos2(1,:)],[],1);
posy = reshape([pos1(2,:);pos2(2,:)],[],1);
posz = reshape([pos1(3,:);pos2(3,:)],[],1);
pos_new = zeros(3*size(pos1,2),3);
ind_NaN = 3:3:size(pos_new,1);
pos_new(ind_NaN,:) = nan(length(ind_NaN),3);
pos_new(~isnan(pos_new(:,1)),:) = [posx posy,posz];
end