function h = bacFigure(y_re,tref,t_str1,leg_str,new_sel)
    if nargin < 5
        new_sel=1;
    end
    if new_sel==1
        h = figure;
    else
        h=[];
    end
    t_str = cellfun(@(x) strjoin(strsplit(x,'-'),''),tref,'UniformOutput',false);
    T = length(t_str);
    
    plot(y_re,'LineWidth',2);
    set(gca,'xlim',[0,T]);
    set(gca,'XTick',floor(linspace(1,T,15)));
    set(gca,'XTickLabel',t_str(floor(linspace(1,T,15))));
    set(gca,'XTickLabelRotation',90)    
    setpixelposition(h,[223,365,1345,420]);
    %legend(leg_str,'Location','best')
    box off
    if ~isempty(t_str1)
        t_str1 = strrep(t_str1,'_','-');
        title(t_str1)
    end
    if ~isempty(leg_str)
        leg_str = cellfun(@(x) strrep(x,'_','-'),leg_str,'UniformOutput',false);
        legend(leg_str,'location','best')
    end
end