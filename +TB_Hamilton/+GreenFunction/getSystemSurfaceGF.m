function gr = getSystemSurfaceGF(syst,TranslationSymmetry,E,k)
arguments
    syst                                    TB_Hamilton.Builder;
    TranslationSymmetry                     (1,:)   double;
    E                                       (1,1)   double;
    k                                       (1,1)   =[];
end
nodes_tag = string(syst.system_graph.Nodes.Name);
nodes_tag_Trans = TB_Hamilton.Tag.get_tag_TranslationSymmetry(nodes_tag,TranslationSymmetry);
H0_tag = setdiff(nodes_tag,nodes_tag_Trans);
V1_tag = TB_Hamilton.Tag.get_tag_TranslationSymmetry(H0_tag,TranslationSymmetry);
ind0 = findnode(syst.system_graph,H0_tag);
ind1 = findnode(syst.system_graph,V1_tag);

if isa(syst.Hamilton,"double")
    H0 = syst.Hamilton(ind0,ind0);  V1 = syst.Hamilton(ind0,ind1);
    gr = TB_Hamilton.Green_Function.SurfaceGreenFunction(H0,V1,E);
elseif isa(syst.Hamilton,"function_handle")
    Ham = syst.Hamilton(k);
    H0 = Ham(ind0,ind0);  V1 = Ham(ind0,ind1);
    gr = full(TB_Hamilton.Green_Function.SurfaceGreenFunction(H0,V1,E));
end

end

