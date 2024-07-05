function tag_TranslationSymmetry = get_tag_TranslationSymmetry(tag,TranslationSymmetry)
    arguments
        tag    (1,:);
        TranslationSymmetry (1,:) ;
    end
    matches = regexp(tag, '\((.*?)\)', 'tokens'); % 使用正则表达式匹配括号中的内容
    if length(matches) ~= 1
        coordinates_str = cellfun(@(x) x{1},matches); % 获取匹配到的内容
    else
        coordinates_str = cellfun(@(x) x,matches);
    end
    coordinates = cellfun(@(x) str2num(strrep(x, ',',' '))'+TranslationSymmetry',...
        coordinates_str,'UniformOutput',false);
    coordinates = round([coordinates{:}]);
    coordinates = cellstr(string(coordinates(1,:))+"," + string(coordinates(2,:)));
    tag_TranslationSymmetry = strrep(reshape(tag,[],1),reshape(coordinates_str,[],1),reshape(coordinates,[],1));
    % tag_TranslationSymmetry is (:,1) vector.
end