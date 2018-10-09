function [ sensors ] = WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,axislimits )

figure(hg); clf;
if strcmp(winddir,'u'); hi = 1; else hi = 2; end; % define wind direction
plot(Y{hi}(:),X{hi}(:),'.','DisplayName',['Grid points']); hold on;
plot(BCs{hi}(:,1),BCs{hi}(:,2),'blackx','DisplayName','BCs'); % plot BCs
for j = 1:length(Wp.turbine.Cry)
    plot([Wp.turbine.Cry(j)+ [-0.5, 0.5]*Wp.turbine.Drotor]',[Wp.turbine.Crx(j)*[1 1]],'DisplayName',['Turbine' num2str(j)],'LineWidth',3);
    text(Wp.turbine.Cry(j),Wp.turbine.Crx(j),num2str(j));
end

if length(sensors{hi}.grid)>0
    plot(y{hi}(sensors{hi}.grid(:,2)),x{hi}(sensors{hi}.grid(:,1)),'ro','DisplayName','Sensors');
end;
lg = legend('-DynamicLegend');
lg.Position(1:2) = [0.77 0.7];
ylabel('x-direction (m)'); xlabel('y-direction (m)');
title(['Grid layout and sensor locations']);
axis equal; 
zoom reset; set(gca,{'xlim','ylim'},axislimits); % restore zoom
drawnow;


% Create pop-up menu for wind direction selection
txt = uicontrol('Style','text',...
    'Position',[30 517 80 20],...
    'String','Wind direction');
if strcmp(winddir,'u')
    popup = uicontrol('Style','popup',...
        'String', {'u','v'},...
        'Position', [110 520 50 20],...
        'Callback', @selwinddir);
else
    popup = uicontrol('Style','popup',...
        'String', {'v','u'},...
        'Position', [110 520 50 20],...
        'Callback', @selwinddir);
end;


% Create push button for adding/removing points
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Add/remove by point',...
    'Position', [30 475 133 20],...
    'Callback', @addremsensors_point);
btn2 = uicontrol('Style', 'pushbutton', 'String', 'Add/remove by range',...
    'Position', [30 450 133 20],...
    'Callback', @addremsensors_range);
btn3 = uicontrol('Style', 'pushbutton', 'String', 'Clear all',...
    'Position', [30 425 133 20],...
    'Callback', @remallsensors);
btn4 = uicontrol('Style', 'pushbutton', 'String', 'Send to workspace',...
    'Position', [30 125 133 20],...
    'Callback', @sendtoworkspace);
btn5 = uicontrol('Style', 'pushbutton', 'String', 'Save file',...
    'Position', [30 100 133 20],...
    'Callback', @savetemplate);
btn6 = uicontrol('Style', 'pushbutton', 'String', 'Load file',...
    'Position', [30 75 133 20],...
    'Callback', @loadtemplate);

    function selwinddir(source,callbackdata)
        winddir = (source.String(source.Value));
        WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,get(gca,{'xlim','ylim'})); % no new sensors generated
    end

    function addremsensors_point(source,callbackdata)
        n = inputdlg('How many points would you like to add/remove? (0 for exit)','Select sensor locations');
        n = str2num(n{1});
        for i = 1:n
            [ym,xm] = ginput(1); % gather input
            [~,js] = min(abs(xm-x{hi}));
            [~,hs] = min(abs(ym-y{hi}));
            
            % add sensor points to the matrix
            if size(sensors{hi}.grid,1) == 0
                j = 0;
            else
                for j = 1:size(sensors{hi}.grid,1)
                    if sensors{hi}.grid(j,:) == [js, hs] % remove sensor points
                        sensors{hi}.grid = [sensors{hi}.grid(1:j-1,:); sensors{hi}.grid(j+1:end,:)];
                        j = -1;
                        break;
                    end;
                end;
            end;
            if j ~= -1
                sensors{hi}.grid = [sensors{hi}.grid; js hs];
            end;
            sensors = WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,get(gca,{'xlim','ylim'}) ); % Update figure
        end;
    end

    function addremsensors_range(source,callbackdata)
        n = inputdlg('How many rectangular areas would you like to select? (0 for exit)','Select sensor locations');
        n = str2num(n{1});
        for i = 1:n
            rect   = getrect;
            jw = min(find((y{hi}-rect(1)>0))); % most western grid point
            js = min(find((x{hi}-rect(2)>0))); % most southern grid point
            je = max(find((rect(1)+rect(3)-y{hi})>0)); % most eastern grid point
            jn = max(find((rect(2)+rect(4)-x{hi})>0)); % most nothern grid point
            
            for r = js:jn
                for s = jw:je
                    if size(sensors{hi}.grid,1) == 0 % add grid points to the matrix
                        j = 0;
                    else
                        for j = 1:size(sensors{hi}.grid,1)
                            if sensors{hi}.grid(j,:) == [r s] % remove sensor points if existant
                                sensors{hi}.grid = [sensors{hi}.grid(1:j-1,:); sensors{hi}.grid(j+1:end,:)];
                                j = -1;
                                break;
                            end;
                        end;
                    end;
                    if j ~= -1
                        sensors{hi}.grid = [sensors{hi}.grid; r s];
                    end;
                end;
            end;
            sensors = WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,get(gca,{'xlim','ylim'}) ); % Update figure
        end;
    end

    function remallsensors(source,callbackdata)
        sensors{hi}.grid = [];
        sensors = WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,get(gca,{'xlim','ylim'}) );
    end

    function postprocessing
        for hi = 1:2
            sensors{hi}.loc=[]; sensors{hi}.obsid=[]; sensors{hi}.N=[];
            for j = 1:size(sensors{hi}.grid,1) % Determine physical sensor locations
                sensors{hi}.loc(j,:) = [x{hi}(sensors{hi}.grid(j,1)), y{hi}(sensors{hi}.grid(j,2))];
            end;
            
            % Convert grid coordinates to observability vector entry numbers
            if hi==1; winddir = 'u'; else; winddir = 'v'; end;
            sensors{hi}.obsid = WFObs_s_sensors_grid2nr( sensors{hi}.grid, Wp, winddir);
            
            % Sort sensor locations numerically
            [sensors{hi}.obsid,b] = sort(sensors{hi}.obsid);
            sensors{hi}.loc       = sensors{hi}.loc(b,:);
            sensors{hi}.grid      = sensors{hi}.grid(b,:);
            sensors{hi}.N         = length(sensors{hi}.obsid);
            sensors{hi}.winddir   = winddir;
        end;
    end
    
    function sendtoworkspace(source,callbackdata)
        if (size(sensors{1}.grid,1) < 1 | size(sensors{2}.grid,1) < 1)
            errorstring = 'You must have at least one sensor for both u and v.';
            errordlg(errorstring); disp(errorstring);
        else
            postprocessing;
            assignin('base','sensors',sensors);
        end
    end

    function savetemplate(source,callbackdata)
        % Post processing
        if (size(sensors{1}.grid,1) < 1 | size(sensors{2}.grid,1) < 1)
            errorstring = 'You must have at least one sensor for both u and v.';
            errordlg(errorstring); disp(errorstring);
        else
            postprocessing; % format output correctly
            Wp_templ = Wp;
            uisave({'sensors','Wp_templ'},['sensors_' Wp.name]);
        end;
    end

    function loadtemplate(source,callbackdata)
        clear sensors
        Wp_templ = [];
        uiopen('*.mat');
        sensorstemp = sensors; clear sensors;
        sensors{1}.grid = sensorstemp{1}.grid;
        sensors{2}.grid = sensorstemp{2}.grid;
        clear sensorstemp;
        if (sum(abs(Wp.mesh.Lx-Wp_templ.mesh.Lx)) > 1e-5 || ...
            sum(abs(Wp.mesh.Ly-Wp_templ.mesh.Ly)) > 1e-5 || ...
            sum(abs(Wp.turbine.Crx-Wp_templ.turbine.Crx)) > 1e-5 || ...
            sum(abs(Wp.turbine.Cry-Wp_templ.turbine.Cry)) > 1e-5)
            disp('WARNING: Meshing from file does not match current mesh.');
        end
        if(sum(abs(Wp.mesh.Nx-Wp_templ.mesh.Nx)) > 1e-5 || ...
           sum(abs(Wp.mesh.Ny-Wp_templ.mesh.Ny)) > 1e-5)
            error('Loaded file is incompatible (Nx/Ny do not match).');
        end;
        WFObs_s_sensors_plot( Wp,BCs,Y,X,y,x,winddir,sensors,hg,get(gca,{'xlim','ylim'}) ); % Update figure
    end
end