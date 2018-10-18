h = gcf;

% Save text and simple objects first
oldCmap = h.Colormap;
colormap(white);
export_fig out.pdf -m6 -transparent -nocrop

% Save bitmap
colormap(oldCmap);
for j = 1:length(h.Children)
    hc = h.Children(j);
    hc.Box = 'off';
    hc.XColor = [1 1 1];
    hc.YColor = [1 1 1];    
    if strcmp(hc.Type,'colorbar')
        hc.Ticks = []; % Remove ticks
    end
    if strcmp(hc.Type,'axes')
        axes(hc);
        hc.XTick = []; % Remove ticks
        hc.YTick = []; % Remove ticks
        hc.XLabel.String = '';
        hc.YLabel.String = '';
        hc.Title.String  = '';
        for ji = 1:length(hc.Children)
            hci = hc.Children(ji);
            if strcmp(hci.Type,'contour') == false
                hci.Visible = 'off';
            end
        end
    end
end
export_fig out.png -m6 -transparent -nocrop