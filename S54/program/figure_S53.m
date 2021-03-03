function h = figure_S53(y_re,tref,t_str1,new_fig)
    if nargin < 4
        new_fig=1;
    end
    t_str = cellfun(@(x) strjoin(strsplit(x,'-'),''),tref,'UniformOutput',false);
    T = length(t_str);
    lim_cut = min(15,T);
    if new_fig>0
        h = figure;
    else
        h=[];
    end
    plot(y_re,'LineWidth',2);
    set(gca,'xlim',[0,T]);
    set(gca,'XTick',floor(linspace(1,T,lim_cut)));
    set(gca,'XTickLabel',t_str(floor(linspace(1,T,lim_cut))));
    set(gca,'XTickLabelRotation',90) 
    if new_fig>0
    setpixelposition(h,[223,365,1345,420]);
    end
    %legend(leg_str,'Location','best')
    box off
    if ~isempty(t_str1)
        title(t_str1)
    end
end