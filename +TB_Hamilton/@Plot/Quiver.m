function Quiver(pos,vec,options)
    arguments
        pos     (:,2)   double;
        vec     (:,2)   double;
        options.?matlab.graphics.chart.primitive.Quiver;
    end
    propertyCell = namedargs2cell(options);
    if isempty(propertyCell)
        propertyCell = {'LineWidth',1};
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
    quiver(pos_(:,1),pos_(:,2),vec(:,1),vec(:,2),propertyCell{:});
end