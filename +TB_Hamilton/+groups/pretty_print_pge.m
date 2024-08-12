function name = pretty_print_pge(g)
% PRETTY_PRINT_PGE returns a human readable string representation of
% PointGroupElement

arguments
    g  TB_Hamilton.groups.PointGroupElement;
end
R = g.R;

if size(R,1)==1
    if R(1,1) == 1
        rot_name = '1';
    else
        rot_name = '-1';
    end
elseif size(R,1) == 2
    if abs(det(R)-1)<1e-6
        theta = atan2(R(2,1),R(1,1));
        if abs(theta) <1e-6
            rot_name = '1';
        else
            rot_name = "R"+"("+name_angle(theta)  +")";
        end
    else
        % mirror
        [V,~] = eig(R);
        n = V(:,1);
        rot_name = "M"+"("+name_axis(n)  +")";
    end

else
    if abs(det(R)-1)<1e-6
        % pure rotation
        [n,theta] = TB_Hamilton.groups.rotation_to_angle(R);

        if abs(theta)<1e-6
            rot_name = "1";
        else
            axis_name = name_axis(n);
            rot_name = "R"+"("+name_angle(theta)+ "," + axis_name +")";
        end
    else
        [n,theta] = TB_Hamilton.groups.rotation_to_angle(-R);

        if abs(abs(theta)-pi)<1e-4
            axis_name = name_axis(n);
            rot_name = "M"+"("+  axis_name +")";
        elseif abs(theta)<1e-4
            rot_name = "I";
        else
            % The S only represents the rotoinversion S(theta,n) = R(theta,n)*I
            axis_name = name_axis(n);
            rot_name = "S"+"("+name_angle(theta/2)+ "," + axis_name +")";
        end
    end
end

if g.conjugate && (~g.antisymmetry)
    ax_name = "T";
elseif g.conjugate && g.antisymmetry
    ax_name = "P";
elseif ~g.conjugate && g.antisymmetry
    ax_name = "C";
else
    ax_name = "";
end
if rot_name == "1" && ax_name ~= ""
    name = ax_name;
else
    name = ax_name + rot_name;
end

if ~isempty(g.U)
    name = {num2str(name) ; ['U=',mat2str(round(real(g.U),2) + round(imag(g.U),2)*1i)]};
end

end

%%%%%%%%% simplify angle
function angle = name_angle(theta)
[num,den] = rat(theta/pi);
if max(den)>100
    [num,den] = rat(round(theta/pi,2));
end
if abs(num) == 1 && den~=1
    angle = char(string(num)+"pi"+"/"+string(den));
    angle(find(angle=='1',1,"first")) = '';
elseif den==1 && num==1
    angle = "pi";
elseif den==1 && num==-1
    angle = "-pi";
else
    angle = string(num)+"*pi"+"/"+string(den);
end
end

%%%%%%%%%%%%%%%%%%%%
function axis_name = name_axis(axis)
axis(abs(axis)<1e-4)=0;
axis1 = axis/abs(min(axis(axis~=0)));
if max(abs(axis1-round(axis1)))<1e-6 && max(abs(axis1))<5
    axis = axis1;
    [num,den] = rat(axis,2);
    axis_name = "[";
    for j1 = 1:length(axis)
        if num(j1) == 0
            axis_name = axis_name + "0" + " ";
        elseif den(j1) == 1
            axis_name = axis_name + string(num(j1)) + " ";
        else
            axis_name = axis_name + string(num(j1)) + "/" + string(den(j1))+" ";
        end
    end
    axis_name = strip(axis_name,'right',' ') +"]";
else
    axis = round(axis,2);
    axis_name = string(mat2str(axis));
end
end