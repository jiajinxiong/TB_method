function T = Density(pos,rho,options)
arguments
    pos     (:,2) double;
    rho     (:,1) double;
    options.?matlab.graphics.chart.primitive.Contour;
end
propertyCell = namedargs2cell(options);
pos_ = unique(pos,"rows");
% Id = [];s
rho1 = zeros(size(pos_,1),1);
if length(pos_) ~= length(pos)
    for j1 = 1:size(pos_,1)
        id = ismember(pos,pos_(j1,:),"rows");
        rho1(j1) = sum(rho(id));
        % Id = [Id,find(id)];
    end
    rho = rho1;
    pos = pos_;
end

x = pos(:, 1);
y = pos(:, 2);
z = rho;

nx = fix(sqrt(length(x)));
F = scatteredInterpolant(x, y, z, 'natural', 'none');
xq = linspace(min(x), max(x), 2*nx);
yq = linspace(min(y), max(y), 2*nx);
[Xq, Yq] = meshgrid(xq, yq);
Zq = F(Xq, Yq);
B = boundary(x, y);
in = inpolygon(Xq, Yq, x(B), y(B));
Xq(~in) = NaN;
Yq(~in) = NaN;
Zq(~in) = NaN;

T = figure("Units","inches","OuterPosition",[.2,.2,10,10],"PaperPosition",[.2,.2,10,10],"Colormap",slanCM(104));
ax = axes("Parent",T,"Position",[.1,.1,.8,.8],'LineWidth',1,'FontSize',13,'FontName','Times New Roman','Box','off');

contourf(ax,Xq, Yq, Zq,256,propertyCell{:});
% colormap(ax,slanCM(104))
colorbar;
end