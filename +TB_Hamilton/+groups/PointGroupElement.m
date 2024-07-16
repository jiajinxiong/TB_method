classdef PointGroupElement
    % POINTGROUPELEMENT is an element of a point group
    %
    % Parameters
    % ------
    % R : Real space rotation action of the operator. Square matrix with
    % size of the number of spatial demension.
    % conjugate : boolean (default false)
    %
    properties
        R;
        conjugate;
        antisymmetry;
        U;
    end

    methods
        function obj = PointGroupElement(R,options)
            arguments
                R           (:,:) double;
                options.conjugate (1,1) logical  = false;
                options.antisymmetry (1,1) logical  = false;
                options.U (:,:)  = [];
            end
            [obj.R,obj.conjugate,obj.antisymmetry,obj.U] = deal(real(R),options.conjugate,options.antisymmetry,options.U);
        end
        function g = mtimes(obj,g2)
            % MTIMES redefines the times
            arguments
                obj TB_Hamilton.groups.PointGroupElement
                g2 TB_Hamilton.groups.PointGroupElement
            end
            g1 = obj;
            [R1,c1,a1,U1] = deal(g1.R,g1.conjugate,g1.antisymmetry,g1.U);
            [R2,c2,a2,U2] = deal(g2.R,g2.conjugate,g2.antisymmetry,g2.U);
            if isempty(U1) | isempty(U2)
                g_U = [];
            elseif c1
                g_U = U1 * conj(U2);
            else
                g_U = U1 * U2;
            end
            g_R = R1 * R2;
            g = TB_Hamilton.groups.PointGroupElement(g_R,conjugate = c1^c2,antisymmetry = a1^a2,U=g_U);
        end
        function g = mpower(obj,num)
            % MPOWER redefine the power
            arguments
                obj TB_Hamilton.groups.PointGroupElement
                num double;
            end
            g = obj.identity();
            if num >=0
                g1 = obj;
            else
                g1 = obj.inv();
            end
            for i_ = 1:abs(num)
                g = g * g1;
            end
        end

        function result = disp(obj)
            strpg = TB_Hamilton.groups.pretty_print_pge(obj);

            sympref('AbbreviateOutput',false);
            sympref('MatrixWithSquareBrackets',true);
            result = symmatrix(char(strpg));
            disp(result);
        end

        function g0 = inv(obj)
            % INV refine the PointGroupElement
            g1 = obj;
            [R1,c1,a1,U1] = deal(g1.R,g1.conjugate,g1.antisymmetry,g1.U);
            if isempty(U1)
                Uinv = [];
            elseif c1
                Uinv = conj(inv(U1));
            else
                Uinv = inv(U1);
            end
            Rinv = inv(R1);
            g0 = TB_Hamilton.groups.PointGroupElement(Rinv,"conjugate",c1,"antisymmetry",a1,"U",Uinv);
        end

        function g0 = identity(obj)
            % IDENTITY return identity element with the same structure as obj.
            dim = length(obj.R);
            R0 = eye(dim);
            if ~isempty(obj.U)
                U0 = eye(size(obj.U));
            else
                U0 = [];
            end
            g0 = TB_Hamilton.groups.PointGroupElement(R0,"conjugate",false,"antisymmetry",false,"U",U0);
        end

    end
end