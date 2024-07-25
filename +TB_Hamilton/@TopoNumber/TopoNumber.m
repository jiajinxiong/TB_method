classdef TopoNumber
    methods(Static)
        function chern = ChernByWilsonLoop(mesh,syst,bandOrEnergy,method)
            arguments
                mesh                pde.FEMesh;
                syst                  TB_Hamilton.Builder;
                bandOrEnergy                 (1,1)   double;
                method   (1,:) char {mustBeMember(method,{'SigleBand','MultiBands','Energy'})} = 'MultiBands';
            end
            points = mesh.Nodes';       ind = mesh.Elements';
            ind_len = length(ind);      chern = 0;
            Hamilton = syst.Hamilton;

            for j1 = 1:ind_len
                kxs = points(ind(j1,:),1);    kys = points(ind(j1,:),2);
                Hn1 = full(Hamilton(kxs(1),kys(1)));
                Hn2 = full(Hamilton(kxs(2),kys(2)));
                Hn3 = full(Hamilton(kxs(3),kys(3)));
                [un1,E1] = eig(Hn1); [un2,E2] = eig(Hn2); [un3,E3] = eig(Hn3);
                if method == "SigleBand"
                    Psi_n = det(un1(:,bandOrEnergy)'*un2(:,bandOrEnergy))*det(un2(:,bandOrEnergy)'*un3(:,bandOrEnergy))*...
                        det(un3(:,bandOrEnergy)'*un1(:,bandOrEnergy));
                elseif method == "MultiBands"
                    Psi_n = det(un1(:,1:bandOrEnergy)'*un2(:,1:bandOrEnergy))*det(un2(:,1:bandOrEnergy)'*un3(:,1:bandOrEnergy))*...
                        det(un3(:,1:bandOrEnergy)'*un1(:,1:bandOrEnergy));
                else
                    energy = bandOrEnergy;
                    num1 = length(find(diag(E1) < energy));
                    num2 = length(find(diag(E2) < energy));
                    num3 = length(find(diag(E3) < energy));
                    num = max([num1,num2,num3]);
                    Psi_n = det(un1(:,1:num)'*un2(:,1:num))*det(un2(:,1:num)'*un3(:,1:num))*...
                        det(un3(:,1:num)'*un1(:,1:num));
                end
                chern = chern -1i*log(Psi_n/abs(Psi_n));
            end
            chern = chern/(2*pi);
        end

        function C = Curvature(syst,kk,num,options)
            arguments
                syst    TB_Hamilton.Builder;
                kk      (2,:)   double;
                num     double;
                options.num_loop        double = 10;
                options.ismulband       double = 0;
            end
            Ham = syst.Hamilton;    nk = size(kk,2);
            Delta = 1e-2;       num_loop = options.num_loop;
            C = zeros(1,nk);
            for j1 = 1:nk
                kx = kk(1,j1);  ky = kk(2,j1);
                Wilson_loop = 1;
                k_loop = [kx+linspace(0,Delta,num_loop),(kx+Delta)*ones(1,num_loop),...
                    kx+linspace(Delta,0,num_loop),kx*ones(1,num_loop);...
                    ky*ones(1,num_loop),ky+linspace(0,Delta,num_loop),...
                    (ky+Delta)*ones(1,num_loop),ky+linspace(Delta,0,num_loop)];

                for j2 = 1:size(k_loop,2)-1
                    H1 = full(Ham(k_loop(1,j2),k_loop(2,j2)));
                    H2 = full(Ham(k_loop(1,j2+1),k_loop(2,j2+1)));
                    [V1,~] = eig(H1);   [V2,~] = eig(H2);
                    if ~options.ismulband
                        Wilson_loop = Wilson_loop * V1(:,num)'*V2(:,num);
                    else
                        Wilson_loop = Wilson_loop * det(V1(:,1:num)'*V2(:,1:num));
                    end
                end
                H1 = full(Ham(k_loop(1,end),k_loop(2,end)));
                H2 = full(Ham(k_loop(1,1),k_loop(2,1)));
                [V1,~] = eig(H1);   [V2,~] = eig(H2);
                if ~options.ismulband
                    Wilson_loop = Wilson_loop * V1(:,num)'*V2(:,num);
                else
                    Wilson_loop = Wilson_loop * det(V1(:,1:num)'*V2(:,1:num));
                end
                C(j1) = 1i*log(Wilson_loop)/Delta^2;
            end
        end

        function [Z_xspin2,Z_yspin1,Es] = Zeeman_curvature(syst,kx,ky,spin1,spin2,BandorEnergy,methods)
            % ZEEMAN_CURVATURE is used to calculate the zeeman berry curvature
            % Z_{x spin2}^n = 2Re sum_m (ir^x_{nm} spin2_{mn}).
            % 
            % Parameters
            % ------
            % syst : TB_Hamilton.Builder.
            % kx : wave vector in x direction.
            % ky : wave vector in y direction.
            % spin1 : spin operator.
            % spin2 : spin operator.
            % BandorEnergy : band index or fermi energy.
            % methods : 'Band' or 'Energy' (default 'Band')
            % 
            % Return
            % ------
            % Z_xspin2 : Zeeman Berry  curvature
            %   defined as Z_{x spin2} = 2Re sum (ir^x_{nm} spin2_{mn})
            %   where `spin2` is the input parameter.
            % Z_yspin1 : Zeeman Berry curvature
            %   defined as Z_{y spin1} = 2Re sum (ir^y_{nm} spin1_{mn})
            %   where `spin1` is the input paramater
            arguments
                syst    TB_Hamilton.Builder;
                kx      (1,:) double;
                ky      (1,:) double;
                spin1   double;
                spin2   double;
                BandorEnergy (1,:) double;
                methods (1,:) char {mustBeMember(methods,{'Band','Energy'})} = 'Band';
            end
            nx = length(kx); ny = length(ky);
            dx = kx(2) - kx(1);
            dy = ky(2) - ky(1);
            dim = syst.system_graph.numnodes;
            if methods == "Band"
                nB = length(BandorEnergy);
                Z_xspin2 = zeros(nx,ny,nB);               Z_yspin1 = zeros(nx,ny,nB);
                EE = zeros(nx,ny,nB);
                parfor j1 = 1:nx
                    for j2 = 1:ny
                        H = full(syst.Hamilton(kx(j1),ky(j2)));
                        [V,E] = eig(H); E = diag(E);

                        H_dx = full(syst.Hamilton(kx(j1)+dx,ky(j2)));
                        H_dy = full(syst.Hamilton(kx(j1),ky(j2)+dy));
                        vx = (H_dx-H)/dx;  vy = (H_dy-H)/dy;
                        for ind = 1:nB
                            nband = BandorEnergy(ind);
                            Z_xy_ = 0; Z_yx_ = 0;
                            for i1 = [1:ind-1,ind+1:dim]
                                Delta = E(nband)-E(i1);
                                Z_xy_ = Z_xy_ + (V(:,nband)'*vx*V(:,i1) * V(:,i1)'*spin2*V(:,nband) )* (dx*dy)/Delta;
                                Z_yx_ =Z_yx_ + (V(:,nband)'*vy*V(:,i1)*V(:,i1)' * spin1 * V(:,nband))* (dx*dy)/Delta;
                            end
                            Z_xspin2(j1,j2,ind) = Z_xy_ + Z_xy_'; Z_yspin1(j1,j2,ind) = Z_yx_ + Z_yx_';
                            EE(j1,j2,ind) = E(nband);
                        end
                    end
                end
            else
                Z_xspin2 = zeros(nx,ny);               Z_yspin1 = zeros(nx,ny);
                parfor j1 = 1:nx
                    for j2 = 1:ny
                        H = full(syst.Hamilton(kx(j1),ky(j2)));
                        [V,E] = eig(H); E = diag(E);
                        H_dx = full(syst.Hamilton(kx(j1)+dx,ky(j2)));
                        H_dy = full(syst.Hamilton(kx(j1),ky(j2)+dy));
                        vx = (H_dx-H)/dx;  vy = (H_dy-H)/dy;
                        nE = length(find(E<=BandorEnergy));
                        Z_xy_ = 0; Z_yx_ = 0;
                        for ind = 1:nE
                            nband = ind;
                            
                            for i1 = [1:ind-1,ind+1:dim]
                                Delta = E(nband)-E(i1);
                                Z_xy_ = Z_xy_ + (V(:,nband)'*vx*V(:,i1) * V(:,i1)'*spin2*V(:,nband) )* (dx*dy)/Delta;
                                Z_yx_ =Z_yx_ + (V(:,nband)'*vy*V(:,i1)*V(:,i1)' * spin1 * V(:,nband))* (dx*dy)/Delta;
                            end
                            Z_xspin2(j1,j2) = Z_xy_ + Z_xy_'; Z_yspin1(j1,j2) = Z_yx_ + Z_yx_';
                        end
                    end
                end
            end
            if nargout == 3
                Es = EE;
            end
        end

        function [Omega_spin,Es] = spin_curvature(syst,Spin,kx,ky,BandorEnergy,methods)
            % SPIN_CURVATURE is used to calculate the spin berry curvature
            % that defined as Omega^s_n =
            % -imag(<n|Spin*vy|m><m|vx|n>/e_{mn}^2), where `Spin` is the
            % second parameter.
            % 
            % Parameters
            % ------
            % syst : TB_Hamilton.Builder
            % Spin : spin operator
            % kx : wave vector x direction
            % ky : wave vector y direction
            % BandorEnergy : band index or fermi energy
            % methods : 'Band' or 'Energy' (default 'Band')
            %
            % Return
            % ------
            % Omega_spin : spin berry curvature
            %   -imag( <n|s*vy|m><m|vx|n> ) where 's' is the Spin.
            % Es : bulk band energy (optional)
            %   when return value number equal 2
            arguments
                syst    TB_Hamilton.Builder;
                Spin    (:,:) double
                kx      (1,:) double;
                ky      (1,:) double;
                BandorEnergy (1,:) double;
                methods (1,:) char {mustBeMember(methods,{'Band','Energy'})} = 'Energy';
            end
            dim = syst.system_graph.numnodes;
            nx = length(kx); ny = length(ky);
            dy = ky(2)-ky(1);  dx = kx(2) - kx(1);
            assert(length(Spin)==dim,'The dimension of Spin must equal Hamiltonian');
            
            if methods == "Band"
                nB = length(BandorEnergy);
                Omega_spin = zeros(nx,ny,nB);
                EE = zeros(nx,ny,nB);
                parfor j1 = 1:nx
                    for j2 = 1:ny
                        H = syst.Hamilton(kx(j1),ky(j2));  [V,E] = eig(full(H)); E = diag(E);
                        H_dx = syst.Hamilton(kx(j1)+dx,ky(j2)); H_dy = syst.Hamilton(kx(j1),ky(j2)+dy);
                        vy = (H_dy-H)/dy;  vx = (H_dx-H)/dx;
                        for i1 = 1:nB
                            Omega_ = 0;
                            nband = BandorEnergy(i1);
                            for i2 = [1:nband-1,nband+1:dim]
                                Delta = E(nband)-E(i2);
                                Omega_ = Omega_ + V(:,nband)'*(vy*Spin+Spin*vy)*V(:,i2)*V(:,i2)'*vx*V(:,nband)/Delta^2*dx*dy;
                            end
                            Omega_spin(j1,j2,i1) = Omega_;
                            EE(j1,j2,i1) = E(nband);
                        end
                    end
                end
            else
                Omega_spin = zeros(nx,ny);
                parfor j1 = 1:nx
                    for j2 = 1:ny
                        H = syst.Hamilton(kx(j1),ky(j2));  [V,E] = eig(full(H)); E = diag(E);
                        H_dx = syst.Hamilton(kx(j1)+dx,ky(j2)); H_dy = syst.Hamilton(kx(j1),ky(j2)+dy);
                        vy = (H_dy-H)/dy;  vx = (H_dx-H)/dx;
                        nE = length(find(E<=BandorEnergy));
                        Omega_ = 0;
                        for i1 = 1:nE
                            nband = i1;
                            for i2 = [1:nband-1,nband+1:dim]
                                Delta = E(nband)-E(i2);
                                Omega_ = Omega_ + V(:,nband)'*(vy*Spin+Spin*vy)*V(:,i2)*V(:,i2)'*vx*V(:,nband)/Delta^2*dx*dy;
                            end
                            Omega_spin(j1,j2) = Omega_;
                        end
                    end
                end
            end
            Omega_spin = -imag(Omega_spin);
            if nargout == 2
                Es = EE;
            end
        end

        function V_ = systFBZ(syst)
            % SYSTFBZ 
            arguments
                syst    TB_Hamilton.Builder;
            end
            a1 = [syst.vec(1,:),0];   a2 = [syst.vec(2,:),0];
            a3 = [0,0,1];
            volume = abs(dot(a1,cross(a2,a3)));
            b1 = 2*pi*cross(a2,a3)/volume;   b2 = 2*pi*cross(a3,a1)/volume;
            b1 = b1(1:2);   b2 = b2(1:2);
            [x,y] = meshgrid(-3:3,-3:3);
            x = x(:);   y = y(:);
            R_Point = double(x * b1 + y * b2);
            DT = delaunayTriangulation(R_Point);
            [V,r] = voronoiDiagram(DT);
            for j1 = 1:length(r)
                V_ = V(r{j1},:);
                if min(V_(:,1))*max(V_(:,1))<0 && min(V_(:,2))*max(V_(:,2))<0 ...
                        &&max(V_(:))<1e3
                    break;
                end
            end
        end

        function mesh = FBZTri(syst,size_grid)
            arguments
                syst        TB_Hamilton.Builder;
                size_grid   (1,1) double=.5;
            end
            FBZ = TB_Hamilton.TopoNumber.systFBZ(syst);
            region=polyshape(FBZ(:,1),FBZ(:,2));
            tr = triangulation(region);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            geometryFromMesh(model,tnodes,telements);
            mesh=generateMesh(model,"Hmax",size_grid,"GeometricOrder","linear");
        end
    end


end