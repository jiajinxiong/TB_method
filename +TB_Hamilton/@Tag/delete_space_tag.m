function tag_end = delete_space_tag(tag)
    matches = regexp(tag, '\((.*?)\)', 'tokens'); % 使用正则表达式匹配括号中的内容
    coordinates_str = cellfun(@(x) x{1},matches); % 获取匹配到的内容
    tag_end = strrep(tag,"("+ coordinates_str + ")",'');
end

