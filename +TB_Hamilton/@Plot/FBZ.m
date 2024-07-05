function FBZ(syst)
arguments
    syst        TB_Hamilton.Builder;
end

% PLOT_FBZ is used to plot FBZ
% OUTPUT
%   V_ : the points of FBZ including the (0,0) point
%
% INSTRUCTIONS
%   plot_FBZ() : plot the FBZ and reciprocal lattice
%   V_ = plot_FBZ() : get the points of FBZ including (0,0).

a1 = double([syst.vec(1,:),0]);   a2 = double([syst.vec(2,:),0]);
a3 = [0,0,1];
volume = abs(dot(a1,cross(a2,a3)));
b1 = 2*pi*cross(a2,a3)/volume;   b2 = 2*pi*cross(a3,a1)/volume;
b1 = b1(1:2);   b2 = b2(1:2);
R_Point = [];
for j1 = -3:3
    for j2 = -3:3
        R_Point = [R_Point;j1*b1+j2*b2];
    end
end
[Vx,Vy] = voronoi(R_Point(:,1),R_Point(:,2));
Vx(abs(Vx)>max(R_Point(:))) = NaN;
Vy(abs(Vy)>max(R_Point(:))) = NaN;

plot(Vx,Vy,'--','color','k','LineWidth',.5);
hold on
% quiver(0,0,b1(1),b1(2),'LineWidth',2,'AutoScale','off');
% quiver(0,0,b2(1),b2(2),'LineWidth',2,'AutoScale','off');
% scatter(R_Point(:,1),R_Point(:,2),'filled','SizeData',50);
% text(2*b1(1)/3,2*b1(2)/3,...
%     '$b_1$','FontSize',20,'Interpreter','latex');
% text(2*b2(1)/3,2*b2(2)/3,...
    % '$b_2$','FontSize',20,'Interpreter','latex');
xlim([-max(abs([b1(1),b2(1)])),max(abs([b1(1),b2(1)]))]);
ylim([-max(abs([b1(2),b2(2)])),max(abs([b1(2),b2(2)]))]);
xticks([]);yticks([]);
% xlabel('$k_x$','Interpreter','latex','FontSize',15);
% label = ylabel('$k_y$','Interpreter','latex','FontSize',15,'Rotation',0);
% set(label, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
axis equal;
end