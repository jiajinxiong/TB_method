function F = Plot_FG(FG,PType)
    arguments
        FG table;
        PType string = "";
    end
    F = figure("Units","inches","OuterPosition",[.2,.2,10,10],"PaperPosition",[.2,.2,10,10]);
    G = FG.uneq_Feymann_Graph;
    topo_eq_num = FG.Number;
    % num = FG.Number;
    nG = length(G);     nF = ceil(sqrt(nG));
    nFx = nF;           nFy = ceil(nG/nFx);
    w = 0.99/nFx; h = 0.99/nFy;
    for j1 = 1:nFx
        for j2 = 1:nFy
            if j1+(j2-1)*nFx<=nG
                axes("Parent",F,"Position",[(j1-1)*w+5e-3,(nFy-j2)*h+5e-3,w,h],'LineWidth',1,'FontSize',13,'FontName','Times New Roman',...
                'Box','off');
                % plot(G{j1+(j2-1)*nF});
                % axis off;
                G1 = G{j1+(j2-1)*nFx};
                pg = plot(G1,"LineStyle",G1.Edges.LineStyle,'linewidth',2,'EdgeColor','k','MarkerSize',7,'NodeColor','k'...
                        ,'NodeFontSize',15,'NodeFontName','Times New Roman','EdgeAlpha',1,'ArrowSize',G1.Edges.Arrow,'Layout','auto');
                legend(sprintf("Topological Eq. Num: %d",topo_eq_num(j1+(j2-1)*nFx)),'FontSize',13,'FontName','Times New Roman');
                if PType == "G"
                    highlight(pg,'Edges',find(G1.Edges.LineStyle == '-'));
                end
            end
        end
    end
end