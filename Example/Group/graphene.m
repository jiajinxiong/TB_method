clc;clear
C3z = TB_Hamilton.groups.rotation(2*pi/3,[0,0,1],U=eye(2));
Mx = TB_Hamilton.groups.rotation(pi,[1,0,0],"inversion",true,U=[0,1;1,0]);
My = TB_Hamilton.groups.rotation(pi,[0,1,0],"inversion",true,U=eye(2));
T = TB_Hamilton.groups.time_reversal(3);
G = TB_Hamilton.groups.generate_group([C3z,Mx,My]);


GCT = TB_Hamilton.groups.Character_Table(G,false);

syms a b c kx ky kz real;

sx = [0,1;1,0]; sy = [0,-1i;1i,0]; sz = [1,0;0,-1];

a1 = [-1,sqrt(3)]/2; a2 = -[1,sqrt(3)]/2; 
f1 = exp(1i*kx) + exp(1i*dot([kx,ky],a1)) + exp(1i*dot([kx,ky],a2));
fun = [0,f1;conj(f1),0];

% fun = ky*sy;
% fun = sx;
% for j1 = 1:6
%     f1 = simplify(TB_Hamilton.groups.Project_Irr(G,j1,fun,[kx,ky,kz]));
%     if max(abs(subs(f1,[kx,ky,kz],rand(1,3))),[],"all")>0
%         j1
%     end
% end
% T = GCT.Equiv_Class;
% simplify(TB_Hamilton.groups.Project_Irr(G,1,fun,[kx,ky,kz])-fun,'Steps',100)
% TB_Hamilton.groups.Irr_decompose(GCT,T(3,:).*T(6,:),false)


% fun = exp(1i*[kx,ky]*[-1,sqrt(3)]'/2)+exp(1i*[kx,ky]*[-1,-sqrt(3)]'/2)+exp(1i*[kx,ky]*[1,0]')*2;
% taylor(fun,[kx,ky],'Order',3)

TB_Hamilton.groups.find_symmetries(fun,G,[kx,ky,kz])
