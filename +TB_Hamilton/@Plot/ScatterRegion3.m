function T = ScatterRegion3(syst,options1,options2)
arguments
    syst        TB_Hamilton.Builder
    options1.linewidth                       (1,1) double = .05;
    options1.colormap_site                        = slanCM(92);
    options1.colormap_hopping                     = slanCM(13);
    options2.size = .1;
end

T = figure("Units","inches","OuterPosition",[1,1,10,10],"PaperPosition",[1,1,10,10]);
ax = axes("Parent",T,'NextPlot','add','Colormap',options1.colormap_site,"Position",[0,0,1,1]);
set(ax,'Visible','off');

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
    Plot_Sites(ax,site_pos1,site_vec(ind,:),syst.site_family(j2),propertyCell{:});
end

if ~isempty(syst.hopping_pos)
    dim = size(syst.hopping_pos,2)/2;
    pos1 = syst.hopping_pos(:,1:dim);
    pos2 = syst.hopping_pos(:,dim+1:end);
    % pos = get_plot_pos(pos1,pos2);
    hold on
    for j1 = 1:size(pos1,1)
        draw3DBar(pos1(j1,:),pos2(j1,:),options1.linewidth,options1.colormap_hopping);
    end
    
    % plot3(pos(:,1),pos(:,2),pos(:,3),'k');
    % hold on
    % if ~isempty(syst.TranslationSymmetry)
    %     V = adjacency(syst.system_graph);
    %     [x,y] = gplot(V,syst.ham_pos(syst.cell_ind,:));
    %     plot3(x,y,'k','LineWidth',1);
    % end
end
hold off



camlight('headlight');
lighting gouraud;
material shiny;
axis equal;

end

function Plot_Sites(ax,site_pos,site_vec,site_family,options)
arguments
    ax
    site_pos    (:,3)   double;
    site_vec       double;
    site_family;
    options.size = .1;
end
% propertyCell = namedargs2cell(options);
[X,Y,Z] = sphere;

for j1 = 1:size(site_pos,1)
    surf(ax,options.size*X+site_pos(j1,1),options.size*Y+site_pos(j1,2),options.size*Z+site_pos(j1,3),"EdgeColor","none");
end
shading interp

% s = scatter3(site_pos(:,1),site_pos(:,2),site_pos(:,3),'filled',propertyCell{:});
% dtt = s.DataTipTemplate;
% site_family_ = dataTipTextRow('Site Family',repmat(site_family,1,size(site_pos,1)));
% site_vec_1 = dataTipTextRow('a_1',site_vec(:,1)');
% site_vec_2 = dataTipTextRow('a_2',site_vec(:,2)');
% dtt.DataTipRows(1) = site_family_;
% dtt.DataTipRows(2) = site_vec_1;
% dtt.DataTipRows(end) = site_vec_2;
% if size(site_vec,2) == 3
%     site_vec_3 = dataTipTextRow('a_3',site_vec(:,3)');
%     dtt.DataTipRows(end+1) = site_vec_3;
% end


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


function draw3DBar(point1, point2, radius,colorMap)
    % point1 和 point2 是三维空间中的两个点
    % radius 是圆柱体的半径

    % 计算两个点之间的向量
    v = point2 - point1;
    L = norm(v); % 计算向量的长度

    % 生成圆柱体
    [X, Y, Z] = cylinder(radius, 20);
    Z = Z * L; % 调整圆柱体的高度

    % 计算旋转矩阵
    v = v / L; % 归一化向量
    if v(3) == 1
        R = eye(3);
    else
        theta = acos(v(3));
        k = [-v(2), v(1), 0];
        k = k / norm(k);
        K = [0, -k(3), k(2); k(3), 0, -k(1); -k(2), k(1), 0];
        R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
    end

    % 旋转和平移圆柱体
    for i = 1:numel(X)
        P = R * [X(i); Y(i); Z(i)];
        X(i) = P(1) + point1(1);
        Y(i) = P(2) + point1(2);
        Z(i) = P(3) + point1(3);
    end
    surf(X, Y, Z,colorMap, 'FaceColor', 'r', 'EdgeColor', 'none');
    % colormap(colorMap);
end
