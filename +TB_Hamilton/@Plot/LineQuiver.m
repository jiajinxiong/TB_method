function LineQuiver(pos,vec,options)
    arguments
        pos     (:,2)   double;
        vec     (:,2)   double;
        options.?matlab.graphics.chart.primitive.Line;
    end
    pos_ = unique(pos,"rows");
    Id = [];
    if length(pos_) ~= length(pos)
        for j1 = 1:size(pos_,1)
            id = ismember(pos,pos_(j1,:),"rows");
            Id = [Id,find(id)];
        end
        vecx = vec(:,1);    vecy = vec(:,2);
        vec = [reshape(sum(vecx(Id)),[],1) reshape(sum(vecy(Id)),[],1)];
    end
    minx = min(pos_(:,1));   maxx = max(pos_(:,1));
    miny = min(pos_(:,2));   maxy = max(pos_(:,2));
    num = 500;
    [X,Y] = meshgrid(linspace(minx,maxx,num),linspace(miny,maxy,num));
    V1 = griddata(pos_(:,1),pos_(:,2),vec(:,1),X,Y,"natural");
    V2 = griddata(pos_(:,1),pos_(:,2),vec(:,2),X,Y,"natural");
    l = streamslice(X,Y,V1,V2,2);
    set(l,options);
end

