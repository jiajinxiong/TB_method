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
            obj.R(abs(obj.R)<1e-3)=0;
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
            if isempty(U1) || isempty(U2)
                g_U = [];
            elseif c1
                g_U = U1 * conj(U2);
            else
                g_U = U1 * U2;
            end
            g_R = R1 * R2;
            g = TB_Hamilton.groups.PointGroupElement(g_R,conjugate = bitxor(c1,c2),antisymmetry = bitxor(a1,a2),U=g_U);
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

        function eq_g = eq(obj,g)
            arguments
                obj TB_Hamilton.groups.PointGroupElement;
                g   TB_Hamilton.groups.PointGroupElement;
            end
            eq_g = arrayfun(@(x) eq_(x,g),obj);
            function eq_AB = eq_(A,B)
                
                R_eq = abs(max(A.R-B.R,[],'all',"ComparisonMethod","abs"))<1e-4;
                basic_eq = R_eq && (A.antisymmetry==B.antisymmetry) && (A.conjugate==B.conjugate);
                if isempty(A.U)&&isempty(B.U)
                    U_eq = 1;
                elseif bitxor(isempty(A.U),isempty(B.U))
                    U_eq = 0;
                else
                    C = A.inv*B;
                    U1 = C.U;
                    if abs(U1(1,1))<1e-4
                        U_eq = 0;
                    else
                        U_eq = abs(max(U1/U1(1,1)-eye(size(U1)),[],'all',"ComparisonMethod","abs"))<1e-4;
                    end
                end
                eq_AB = all([U_eq, basic_eq]);
            end
        end
        function result = lt(obj,g)
            result = arrayfun(@(x) lt_(x,g),obj);
            function result_ = lt_(A,B)
                if ~isequal([A.conjugate,A.antisymmetry],[B.conjugate,B.antisymmetry])
                    if A.conjugate<B.conjugate
                        result_ = true;
                    elseif (A.antisymmetry < B.antisymmetry) && A.conjugate==B.conjugate
                        result_ = true;
                    else
                        result_ = false;
                    end
                else
                    result_ = string(char(reshape(A.R,1,[])))<string(char(reshape(B.R,1,[])));
                end
            end
        end



        function GT = disp(obj)
            arguments
                obj TB_Hamilton.groups.PointGroupElement;
            end
            strpg = arrayfun(@(x) TB_Hamilton.groups.pretty_print_pge(x,false),obj,'UniformOutput',false);
            strpg = [strpg{:}]';
            GT = table;
            GT.R = strpg(:,1);
            GT.U = strpg(:,2);
        end

        function result = apply(obj,model,k)
            arguments
                obj     TB_Hamilton.groups.PointGroupElement
                model   ;
                k       = [];
            end  
            if ~isempty(k)
                if obj.conjugate
                    result = subs(model,k, -k * obj.R);
                else
                    result = subs(model,k, k * obj.R);
                end
            end
            if obj.conjugate
                result = conj(result);
            end
            if obj.antisymmetry
                result = -result;
            end
            if ~isempty(obj.U)
                result = obj.U * result * obj.U';
            end
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