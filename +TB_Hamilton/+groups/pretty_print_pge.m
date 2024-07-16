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
        axis_name = name_axis(n);
        if abs(theta-pi)<1e-4
            rot_name = "M"+"("+  axis_name +")";
        else
            rot_name = "I*R"+"("+name_angle(theta)+ "," + axis_name +")";
        end
    end
    
end
name = rot_name;
end


function angle = name_angle(theta)
[num,den] = rat(theta/pi);
if abs(num) == 1 && den~=1
    angle = char(string(num)+"pi"+"/"+string(den));

    angle(find(angle=='1',1,"first")) = '';
elseif den == 1 && num~=1
    angle = string(num)+"*pi";
elseif den==1 && num==1
    angle = "pi";
else
    angle = string(num)+"*pi"+"/"+string(den);
end
end

function axis_name = name_axis(axis)
axis1 = axis/abs(min(axis(axis~=0)));
if max(abs(axis1-round(axis1)))<1e-6
    axis = axis1;
end

[num,den] = rat(round(axis,2));
axis_name = "[";
for j1 = 1:length(axis)
    if num(j1) == 0
        axis_name = axis_name + "0" + ",";
    elseif den(j1) == 1
        axis_name = axis_name + string(num(j1)) + ",";
    else
        axis_name = axis_name + string(num(j1)) + "/" + string(den(j1))+",";
    end
end

axis_name = strip(axis_name,'right',',') +"]";
end