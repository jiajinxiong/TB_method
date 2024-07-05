classdef Tag    
    methods(Static)
        deleta_space_tag();
        tag_TranslationSymmetry = get_tag_TranslationSymmetry(tag,TranslationSymmetry);

        function H0_tag = get_H0_tag(syst)
            arguments 
                syst    TB_Hamilton.Builder
            end
            if size(syst.TranslationSymmetry,1) ~= 1
                error("The get_H0_tag can only be applied to semi-infinite system");
            end
            site_tag = string(table2array(syst.system_graph.Nodes));
            Translation_site_tag = TB_Hamilton.Tag.get_tag_TranslationSymmetry(site_tag,syst.TranslationSymmetry);
            H0_tag = setdiff(site_tag,Translation_site_tag);
        end
        function coordinates = get_vec_tag(tag)
            matches = regexp(tag, '\((.*?)\)', 'tokens'); % 使用正则表达式匹配括号中的内容
            coordinates_str = cellfun(@(x) x{1},matches); % 获取匹配到的内容
            coordinates = cellfun(@(x) str2num(strrep(x, ',',' '))',coordinates_str,'UniformOutput',false);
            coordinates = [coordinates{:}];
        end
    end
end

