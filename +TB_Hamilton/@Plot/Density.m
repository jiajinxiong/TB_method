function Density(pos,rho,alpha)
arguments
    pos     (:,2) double;
    rho     (:,1) double;
    alpha   (1,1) double=2;
end
pos_ = unique(pos,"rows");
Id = [];
if length(pos_) ~= length(pos)
    for j1 = 1:size(pos_,1)
        id = ismember(pos,pos_(j1,:),"rows");
        Id = [Id,find(id)];
    end
    rho = reshape(sum(rho(Id)),[],1);
    pos = pos_;
end

shp = alphaShape(pos);
shp.Alpha = alpha;
tri = alphaTriangulation(shp);
% plot(shp);
% figure()
patch('faces',tri,'vertices',pos,'FaceVertexCData',rho,...
    'facecolor','interp','edgecolor','none')
%,'FaceVertexAlphaData',rho ,'FaceAlpha','interp'
axis equal
end