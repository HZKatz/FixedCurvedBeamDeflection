%% ----- SCRIPT INFORMAITON -----
% PURPOSE
% Authors: Harrison Katz - hzkatz@gmail.com
%          Frederick Wachter - wachterfreddy@gmail.com
% Created: 2017-03-25

% Please refer to the Wiki for instructions on how to use this script
% GITHUB WIKI: https://github.com/SOMEPATH

classdef curvedBeam < handle
    
    properties(SetAccess = protected)
        
    end
    
    properties(SetAccess = private, GetAccess = private)
        display
        funct
        param
        react
    end

%% CONSTRUCTOR METHOD
    methods
        
        function solver = curvedBeam(varargin)
        % Construct the structure
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-25
        
            % Parameter Properties
            solver.param.load.x = 0;
            solver.param.load.y = 0;
            solver.param.load.z = 0;
            solver.param.load.alpha = 0;
            solver.param.load.beta  = 0;
            
            solver.param.T = 0;
            solver.param.d = 0;
            solver.param.E = 0;
            solver.param.R = 0;
            solver.param.v = 0;
            
            solver.param.units.force  = [];
            solver.param.units.length = [];
            
            solver.param.I = 0;
            solver.param.J = 0;
            solver.param.G = 0;
            
            solver.param.set = zeros(1,12);
            
            % Function Properties
            solver.funct.I = @(d)(pi*(d^4)/64);
            solver.funct.J = @(d)(pi*(d^4)/32);
            solver.funct.G = @(E,v)(E/(2*(1 + v)));
            solver.funct.load.x = [];
            solver.funct.load.y = [];
            solver.funct.load.z = [];
            
            % Reactions Properties
            solver.react.Mt = 0;
            solver.react.Mr = 0;
            solver.react.Mz = 0;
            solver.react.Fx = 0;
            solver.react.Fy = 0;
            solver.react.Fz = 0;
            
            % Display Properties
            solver.display.figureHandle = [];
            solver.display.axesHandle   = [];
            solver.display.position     = [];
            solver.display.outputIndex  = 1;
            
            solver.initialize();
            
        end
        
    end
    
%% PUBLIC METHODS

    methods
        
        function initialize(solver)
        % Initialize the figure
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-25    
            
            solver.createFigure();
            
        end
        
    end

%% PRIVATE METHODS

    methods (Access = private)
        
        function createFigure(solver)
        % Initialize the figure
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-25
        
            % Create Figure
            solver.display.figureHandle = figure('NumberTitle', 'off', 'Name', 'Curved Beam');
            solver.display.position = get(solver.display.figureHandle, 'Position');
            solver.display.position(4) = solver.display.position(3);
            set(solver.display.figureHandle, 'Position', solver.display.position);
            
            % Display Axes
            solver.display.axesHandle = axes('Units', 'Pixels', 'Position', [50, 210, solver.display.position(3)-90, solver.display.position(4)-230]); plot3(1,1,1);
            xlabel('X-Deformation'); ylabel('Y-Deformation'); zlabel('Z-Deformation');
            title('Curved Beam Deflection'); 
            hold on;

            % Display Control Inputs
            yLoadPosition = 130;
            xFigureMiddle = solver.display.position(3)/2;
            uicontrol('Style', 'text', 'Position', [xFigureMiddle - 175, yLoadPosition + 30, 350, 25], 'String', '________________ Functionality ________________', 'FontSize', 12);
            uicontrol('Style', 'text', 'Position', [xFigureMiddle - 135, yLoadPosition - 5, 40, 25], 'String', 'LOAD', 'FontSize', 12);
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 90, yLoadPosition, 40, 25], 'String', 'X', 'FontSize', 12, 'Callback', {@loadCallback, 1});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 40, yLoadPosition, 40, 25], 'String', 'Y', 'FontSize', 12, 'Callback', {@loadCallback, 2});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 10, yLoadPosition, 40, 25], 'String', 'Z', 'FontSize', 12, 'Callback', {@loadCallback, 3});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 60, yLoadPosition, 40, 25], 'String', 'Alpha', 'FontSize', 12, 'Callback', {@loadCallback, 4});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 110, yLoadPosition, 40, 25], 'String', 'Beta', 'FontSize', 12, 'Callback', {@loadCallback, 5});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 90, yLoadPosition - 32, 40, 25], 'String', 'T', 'FontSize', 12, 'Callback', {@loadCallback, 6});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 40, yLoadPosition - 32, 40, 25], 'String', 'd', 'FontSize', 12, 'Callback', {@loadCallback, 7});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 10, yLoadPosition - 32, 40, 25], 'String', 'E', 'FontSize', 12, 'Callback', {@loadCallback, 8});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 60, yLoadPosition - 32, 40, 25], 'String', 'R', 'FontSize', 12, 'Callback', {@loadCallback, 9});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle + 110, yLoadPosition - 32, 40, 25], 'String', 'v', 'FontSize', 12, 'Callback', {@loadCallback, 10});
            uicontrol('Style', 'text', 'Position', [xFigureMiddle - 139, yLoadPosition - 69, 40, 25], 'String', 'UNITS', 'FontSize', 12);
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 90, yLoadPosition - 64, 65, 25], 'String', 'Force', 'FontSize', 12, 'Callback', {@unitCallback, 1});
            uicontrol('Style', 'edit', 'Position', [xFigureMiddle - 15, yLoadPosition - 64, 65, 25], 'String', 'Length', 'FontSize', 12, 'Callback', {@unitCallback, 2});
            uicontrol('Style', 'text', 'Position', [xFigureMiddle - 155, yLoadPosition - 98, 60, 25], 'String', 'OUTPUT', 'FontSize', 12);
            uicontrol('Style', 'popup', 'Position', [xFigureMiddle - 95, yLoadPosition - 95, 150, 25], 'String', {'Total Deformation', 'X Deformation', 'Y Deformation', 'Z Deformation'}, 'FontSize', 12, 'Callback', @outputCallback);
            uicontrol('Style', 'pushbutton', 'Position', [xFigureMiddle - 152, yLoadPosition - 125, 93, 25], 'String', 'Update', 'FontSize', 12, 'Callback', {@buttonCallback, 1});
            uicontrol('Style', 'pushbutton', 'Position', [xFigureMiddle - 45, yLoadPosition - 125, 93, 25], 'String', 'Calculate', 'FontSize', 12, 'Callback', {@buttonCallback, 2});
            uicontrol('Style', 'pushbutton', 'Position', [xFigureMiddle + 61, yLoadPosition - 125, 93, 25], 'String', 'Animate', 'FontSize', 12, 'Callback', {@buttonCallback, 3});
            
            function loadCallback(src, ~, index)
            % Callback function for text input boxes
            % Author:  Frederick Wachter - wachtefreddy@gmail.com
            % Created: 2017-03-25
            
                value = get(src, 'String');
                if ~(isstrprop(value, 'digit'))
                    solver.param.set(index) = 0;
                    fprintf('[ERROR] Input for text editor box was not an numeric value');
                    % Add error box display message
                    return;
                else
                    value = str2double(value);
                end
                
                switch index
                    case 1 % if x value
                        solver.param.load.x = value;
                    case 2 % if y value
                        solver.param.load.y = value;
                    case 3 % if z value
                        solver.param.load.z = value;
                    case 4 % if alpha value
                        solver.param.load.alpha = value;
                    case 5 % if beta value
                        solver.param.load.beta = value;
                    case 6 % if T value
                        solver.param.T = (value - 90)*(pi/180);
                    case 7 % if d value
                        solver.param.d = value;
                    case 8 % if E value
                        solver.param.E = value;
                    case 9 % if R value
                        solver.param.R = value;
                    case 10 % if v value
                        solver.param.v = value;
                end
                
                solver.param.set(index) = 1;
                
            end
            
            function unitCallback(src, ~, index)
            % Callback function edit boxed used to define units
            % Author:  Frederick Wachter - wachtefreddy@gmail.com
            % Created: 2017-03-27
                
                value = get(src, 'String');
            
                switch index
                    case 1 % if force units definition
                        solver.param.units.force = value;
                        solver.param.set(11) = 1;
                    case 2 % if length units definition
                        solver.param.units.length = value;
                        solver.param.set(12) = 1;
                end
                
            end
            
            function outputCallback(src, ~)
            % Callback function for dropdown menu
            % Author:  Frederick Wachter - wachtefreddy@gmail.com
            % Created: 2017-03-26
            
                solver.display.outputIndex = get(src, 'Value');
                
            end
            
            function buttonCallback(~, ~, index)
            % Callback function for buttons
            % Author:  Frederick Wachter - wachtefreddy@gmail.com
            % Created: 2017-03-26
            
                switch index
                    case 1 % if update button was pressed
                        if (sum(solver.param.set) == 12)
                            cla; fprintf('[INFO] Updating plot\n');
                            solver.loadEquations();
                            solver.loadReactions();
                            solver.getMaterialProperties();
                            solver.update();
                        else
                            fprintf('[ERROR] Not all load parameters have been set');
                            % Add error box display message
                            return;
                        end
                    case 2 % if calculate button was pressed
                        % Perform calculations
                    case 3 % if animate button was pressed
                        % Animate graph
                end
                
            end
            
        end
        
        function update(solver)
        % Update plot
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-27
            
            forceMag = sqrt(solver.param.load.x^2 + solver.param.load.y^2 + solver.param.load.z^2);
            xMag = solver.param.load.x/forceMag;
            yMag = solver.param.load.y/forceMag;
            zMag = solver.param.load.z/forceMag;
            
            meshSize = 51;
            percentAngle = (solver.param.load.alpha + pi/2)/(solver.param.T + pi/2);
            mesh1 = round(meshSize*percentAngle);
            mesh2 = meshSize - mesh1;

            angleArray1 = linspace(-pi/2,solver.param.load.alpha,mesh1);
            angleArray2 = linspace(solver.param.load.alpha,solver.param.T,mesh2);
            angleArray = [angleArray1, angleArray2];

            % Substitute Symoblic Reactions for the solved reactions from user input
            syms Fx Fy Fz Mt Mr Mz
            deflectionX1 = subs(solver.funct.load.x.defQx{1},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);
            deflectionX2 = subs(solver.funct.load.x.defQx{2},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);
            deflectionY1 = subs(solver.funct.load.y.defQy{1},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);
            deflectionY2 = subs(solver.funct.load.y.defQy{2},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);
            deflectionZ1 = subs(solver.funct.load.z.defQz{1},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);
            deflectionZ2 = subs(solver.funct.load.z.defQz{2},[Fx,Fy,Fz,Mt,Mr,Mz],[solver.react.Fx,solver.react.Fy,solver.react.Fz,solver.react.Mt,solver.react.Mr,solver.react.Mz]);

            syms Px Py Pz E I G J R a;
            deflectionX1 = subs(deflectionX1,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            deflectionX2 = subs(deflectionX2,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            deflectionY1 = subs(deflectionY1,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            deflectionY2 = subs(deflectionY2,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            deflectionZ1 = subs(deflectionZ1,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            deflectionZ2 = subs(deflectionZ2,[Px,Py,Pz,E,I,G,J,R,a],[solver.param.load.x,solver.param.load.y,solver.param.load.z,solver.param.E,solver.param.I,solver.param.G,solver.param.J,solver.param.R,solver.param.load.alpha]);
            
            % Define X,Y,Z Coordinates for the undeformed, original beam
            normXCoords = solver.param.R*sin(angleArray);
            normYCoords = solver.param.R*cos(angleArray);
            normZCoords = zeros(1,length(angleArray));

            % Evaluate the deflection of the beam at any angle Beta
            % This is done by substituting the angle array for the variable Beta
            syms b;
            fprintf('[INFO] Performing angle substitution 1/6\n');
            defXCoords1 = eval(subs(deflectionX1, b, angleArray1));
            fprintf('[INFO] Performing angle substitution 2/6\n');
            defXCoords2 = eval(subs(deflectionX2, b, angleArray2));
            fprintf('[INFO] Performing angle substitution 3/6\n');
            defYCoords1 = eval(subs(deflectionY1, b, angleArray1));
            fprintf('[INFO] Performing angle substitution 4/6\n');
            defYCoords2 = eval(subs(deflectionY2, b, angleArray2));
            fprintf('[INFO] Performing angle substitution 5/6\n');
            defZCoords1 = eval(subs(deflectionZ1, b, angleArray1));
            fprintf('[INFO] Performing angle substitution 6/6\n');
            defZCoords2 = eval(subs(deflectionZ2, b, angleArray2));

            defXCoords = [defXCoords1, defXCoords2];
            defYCoords = [defYCoords1, defYCoords2];
            defZCoords = [defZCoords1, defZCoords2];

            % Calculate the total deflection at every point using Pythagorian Thm.
            deflectionMag = sqrt((defXCoords).^2 + (defYCoords).^2 + (defZCoords).^2);

            xScale = 10E1;
            yScale = 10E1;
            zScale = 10E1;

            % xScale = 1; yScale = 1; zScale = 1;

            % Deformed Result is created by adding the deformed coordinates to original
            % The deformed coordinates are scalled to give a better look at def. shape
            xCoords = normXCoords + xScale*defXCoords;
            yCoords = normYCoords + yScale*defYCoords;
            zCoords = normZCoords + zScale*defZCoords;

            % Define a colormapping scale
            switch solver.display.outputIndex
                case 1
                    c = deflectionMag;
                case 2
                    c = defXCoords;
                case 3
                    c = defYCoords;
                case 4
                    c = defZCoords;
            end

            % Create a colormap and save it as a variable
            cmap = colormap;

            % C is now used to index through the colormapping values (min --> max)
            cOut = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
            cOut(isnan(cOut)) = 0;

            % For loop to create line segments with the color map based on deflection
            plot3(normXCoords,normYCoords,normZCoords,'linewidth',4,'color',[.8,.8,.8]);
            
            % Create scatter plots for points used for undeformed and deformed result
            scatter3(normXCoords,normYCoords,normZCoords,7,'k');
            scatter3(xCoords,yCoords,zCoords,7,'k');

            for j = 1:(meshSize-1)
                line(xCoords(j:j+1),yCoords(j:j+1),zCoords(j:j+1),'linewidth',3,'color',cmap(cOut(j),:))
            end

            loadLoc = [xCoords(mesh1), yCoords(mesh1), zCoords(mesh1)];
            scatter3(loadLoc(1), loadLoc(2), loadLoc(3), 8, 'r');
            quiver3(loadLoc(1), loadLoc(2), loadLoc(3), xMag, yMag, zMag, 1, 'linewidth', 2, 'color','r');

            % Define a colorbar and set the limits as min/max deflection magnitude
            h = colorbar;
            caxis([min(c),max(c)]);
            set(get(h,'title'),'string',sprintf('%s',solver.param.units.length));

            % Add labels and change the aspect ratio to better view vertical deflection
            xlabel('X-Deformation'); ylabel('Y-Deformation'); zlabel('Z-Deformation');
            title('Curved Beam Deflection');
            daspect([1,1,1]);
            
            fprintf('\n');
            
        end
        
        function loadEquations(solver)
        % Load the energy method equations
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-27
            
            fprintf('[INFO] Loading energy equations\n');
            solver.funct.load.x = load('deflectX-180deg.mat');
            solver.funct.load.y = load('deflectY-180deg.mat');
            solver.funct.load.z = load('deflectZ-180deg.mat');
            
        end
        
        function loadReactions(solver)
        % Load the reaction equations
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-27
            
            fprintf('[INFO] Loading reaction equations\n');
            reactions = load('reactionForces-Gendeg.mat');
            syms Px Py Pz T a E I G J R;
            
            fprintf('[INFO] Solving reaction equation 1/6\n');
            solver.react.Fx = round(eval(subs(reactions.X.Fx)*100))/100;
            fprintf('[INFO] Solving reaction equation 2/6\n');
            solver.react.Fy = round(eval(subs(reactions.X.Fy)*100))/100;
            fprintf('[INFO] Solving reaction equation 3/6\n');
            solver.react.Fz = round(eval(subs(reactions.X.Fz)*100))/100;
            fprintf('[INFO] Solving reaction equation 4/6\n');
            solver.react.Mt = round(eval(subs(reactions.X.Mt)*100))/100;
            fprintf('[INFO] Solving reaction equation 5/6\n');
            solver.react.Mr = round(eval(subs(reactions.X.Mr)*100))/100;
            fprintf('[INFO] Solving reaction equation 6/6\n');
            solver.react.Mz = round(eval(subs(reactions.X.Mz)*100))/100;
            
        end
        
        function getMaterialProperties(solver)
        % Calculated remaining material properties
        % Author:  Frederick Wachter - wachtefreddy@gmail.com
        % Created: 2017-03-27
        
            fprintf('[INFO] Solving for beam properties\n');
            solver.param.I = solver.funct.I(solver.param.d);
            solver.param.J = solver.funct.J(solver.param.d);
            solver.param.G = solver.funct.G(solver.param.E, solver.param.v);
            
        end
        
    end
    
end


