classdef Monatomic %< TB_Hamilton.Polyatomic
    %MONATOMIC 此处显示有关此类的摘要
    %   此处显示详细说明
    properties
        prim_vecs;
        names;
        basises;
        orbits;
        shape;
    end
    methods
        function obj = Monatomic(prim_vecs, name, offset,orbits)
            arguments
                prim_vecs       double
                name            string
                offset  (1,:)   double;
                orbits  (1,1)   double = 1;
            end
            obj.prim_vecs = prim_vecs;
            obj.names = name;
            obj.basises = offset;
            obj.orbits = orbits;
        end
        function result = subsref(obj, S)
            if strcmp(S(1).type, '()') && max(size(obj)) == 1 ...
                    && isa(obj,'TB_Hamilton.Monatomic')
                result = obj.site_tag(S.subs{:});
            else
                if size({S.subs},2)>1
                    prop = S(1).subs;
                    result = obj.(prop)(S(2).subs{:});
                else
                    prop = S(1).subs;
                    result = obj.(prop);
                end
            end
        end
        function result = site_tag(obj,x,y,z)
            arguments
                obj     TB_Hamilton.Monatomic;
                x       ;
                y       ;
                z       =[];
            end
            if isempty(z)
                if length(x)~=length(y)
                    [x,y] = meshgrid(x,y);
                end
                input = [x(:) y(:)];
            else
                input = [x(:) y(:) z(:)];
            end
                
                len = length(x);
                poses = obj.pos(input);
                str1 = obj.names+"("+join(string(input),',')+")";
                result.pos = poses(1,:);
                result.site_family = obj.names;
                if obj.orbits == 1
                    result.tag = str1(1);
                    for j1 = 2:len
                        result(j1).pos = poses(j1,:);
                        result(j1).tag = str1(j1);
                        result(j1).site_family = obj.names;
                    end
                else
                    result.tag = str1(1)+"_"+ string(1:obj.orbits);
                    for j1 = 2:len
                        result(j1).pos = poses(j1,:);
                        result(j1).tag = str1(j1)+"_"+ string(1:obj.orbits);
                        result(j1).site_family = obj.names;
                    end
                end
                result = result';
        end
        function result = pos(obj,tag)
            arguments
                obj TB_Hamilton.Monatomic
                tag  double
            end
            result = tag * obj.prim_vecs + obj.basises;
        end

        function lat_tag = LatInShape(obj)
            arguments
                obj     TB_Hamilton.Monatomic;
            end
            if ~isa(obj.shape,'polyshape')
                error("Only LatInShape can be used when the shape is well defined.");
            end
            [bx,by] = boundary(obj.shape);
            mn = ceil(max(abs(([bx by] - obj.basises)/obj.prim_vecs)));
            [M,N] = meshgrid(-mn(1):mn(1),-mn(2):mn(2));
            Pointx = M * obj.prim_vecs(1,1) + N * obj.prim_vecs(2,1) + obj.basises(1);
            Pointy = M * obj.prim_vecs(1,2) + N * obj.prim_vecs(2,2) + obj.basises(2);
            in = inpolygon(Pointx,Pointy,bx,by);
            lat_tag = [M(in),N(in)];
        end
    end
end