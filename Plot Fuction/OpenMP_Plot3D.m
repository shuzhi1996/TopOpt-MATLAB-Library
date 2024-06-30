classdef OpenMP_Plot3D
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
% class for Model
% How to use:
%   Type OpenMP_Plot3D.[function name] on MATLAB command window or MATLAB script. 
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

methods (Static)
    %% 单纯的Plot
    function Plot(eta,angle1,angle2,varargin)
        % eta 画图的阈值
        % angle1 和 angle2 旋转的角度（规则详见MATLAB Plot Rotate的官网网页）
        % varargin
        nPlot=max(0,nargin-3);
        for nn = 1:nPlot
            xPlot=varargin;
            
            Aplot2 = xPlot{nn};
            [nely,nelx,nelz] = size(Aplot2);
            ngrid = 1;
            nelyi = ngrid*nely; nelxi = ngrid*nelx; nelzi = ngrid*nelz;
            elyi = nely/(nelyi-1); elxi = nelx/(nelxi-1); elzi = nelz/(nelzi-1);
            ely0 = nely/(nely-1); elx0 = nelx/(nelx-1);   elz0 = nelz/(nelz-1);

            [x,y,z] = ndgrid(0:ely0:nely,0:elx0:nelx,0:elz0:nelz);
            F = griddedInterpolant(x,y,z,Aplot2);
            [X,Y,Z] = ndgrid(0:elyi:nelyi,0:elxi:nelxi,0:elzi:nelzi);
            Aplot2 = F(X,Y,Z);

            Aplot2(:,:,1) = 0;  Aplot2(:,:,end) = 0;
            Aplot2(:,1,:) = 0;  Aplot2(:,end,:) = 0;
            Aplot2(1,:,:) = 0;  Aplot2(end,:,:) = 0;

            Aplot2(1:1,1:1,1:1)=1;
            Aplot2(end,end,end)=1;
            Bplot = Aplot2; permute(Aplot2,[1 3 2]);

            isovals = Bplot; shiftdim(Bplot, 2);
            isovals = smooth3(isovals,'box',1);

            % subplot(Figy,Figx,nn);
            [F1,V1] = isosurface(isovals,eta);
            [F2,V2] = isocaps(isovals,eta);
            F3 = [F1;F2+size(V1,1)];
            V3 = [V1;V2];

            % Sa = isosurface(isovals,eta);
            Sa.Vertices = V3;
            Sa.Faces = F3;
            Sa.FaceColor = [0 0 1];
            Sa.EdgeColor = 'none';
            patch(Sa);
            % patch(isosurface(isovals,eta),'FaceColor',[0 0 1],'EdgeColor','none');
            % patch(isocaps(isovals,eta),'FaceColor',[1 0 0],'EdgeColor','none');
            view([angle1,angle2]); axis equal tight off; camlight; drawnow
        end
    end

    function PlotNode(eta,angle1,angle2,varargin)
        xPlot=varargin;
        Aplot2 = xPlot{1};
        [nely,nelx,nelz] = size(Aplot2);
        [x,y,z] = ndgrid(1:nely,1:nelx,1:nelz);
        F = griddedInterpolant(x,y,z,Aplot2);
        [X,Y,Z] = ndgrid(1:nely+1,1:nelx+1,1:nelz+1);
        Aplot2 = F(X,Y,Z);

        Aplot2(:,:,1) = 0;  Aplot2(:,:,end) = 0;
        Aplot2(:,1,:) = 0;  Aplot2(:,end,:) = 0;
        Aplot2(1,:,:) = 0;  Aplot2(end,:,:) = 0;

        Aplot2(1:1,1:1,1:1)=1;
        Aplot2(end,end,end)=1;
        Bplot = permute(Aplot2,[1 3 2]);

        isovals = shiftdim(Bplot, 2);
        isovals = smooth3(isovals,'box',1);

        [F1,V1] = isosurface(isovals,eta);
        [F2,V2] = isocaps(isovals,eta);
        F3 = [F1;F2+size(V1,1)];
        V3 = [V1;V2];

        % Sa = isosurface(isovals,eta);
        Sa.Vertices = V3;
        Sa.Faces = F3;

        SS = patch(Sa);
        Bplot22 = permute(xPlot{2},[1 3 2]);
        Bplot22 = shiftdim(Bplot22, 2);
        isocolors(Bplot22,SS);
        SS.FaceColor = 'interp';
        SS.EdgeColor = 'none';

        view([angle1,angle2]); axis equal tight off; camlight; drawnow
        colormap(xPlot{3}); colorbar('eastoutside');
    end

    function PlotEle(eta,angle1,angle2,FlagTitle,varargin)
        xPlot=varargin;
        Aplot2 = xPlot{1};

        isovals = Aplot2;
        isovals = smooth3(isovals,'box',1);

        [F1,V1] = isosurface(isovals,eta);
        [F2,V2] = isocaps(isovals,eta);
        F3 = [F1;F2+size(V1,1)];
        V3 = [V1;V2];
        Sa.Vertices = V3;
        Sa.Faces = F3;
        SS = patch(Sa);

        Bplot22 = xPlot{2};
        isocolors(Bplot22,SS);
        SS.FaceColor = 'interp';
        SS.EdgeColor = 'none';

        view([angle1,angle2]); axis equal tight off; camlight; drawnow
        colormap(xPlot{3}); 
%         colorbar('eastoutside');
        ax = gca;
        ax.Interactions = [rotateInteraction dataTipInteraction];

        if FlagTitle > 0
            title(strcat('Itr:  ',num2str(varargin{4}),'  Obj:  ',num2str(round(varargin{5},3)),...
                '  Cons:  ',num2str(round(varargin{6},3)),'  Volfrac:  ',num2str(round(varargin{7},3))))
        end
    end

    %% 翻转功能
    function Aplot2 = flipFilter(Aplot,direction)
        % xz
        if direction == "ymin"
            Aplot2 = zeros(size(Aplot,1)*2,size(Aplot,2)*1,size(Aplot,3));
            Aplot2(1:size(Aplot,1),:,:) = flipud(Aplot);
            Aplot2(size(Aplot,1)+1:end,:,:) = Aplot;
        elseif direction == "ymax"
            Aplot2 = zeros(size(Aplot,1)*2,size(Aplot,2)*1,size(Aplot,3));
            Aplot2(1:size(Aplot,1),:,:) = Aplot;
            Aplot2(size(Aplot,1)+1:end,:,:) = flipud(Aplot);
        elseif direction == "xmin"
            Aplot2 = zeros(size(Aplot,1)*1,size(Aplot,2)*2,size(Aplot,3));
            Aplot2(:,1:size(Aplot,2),:) = flip(Aplot,2);
            Aplot2(:,size(Aplot,2)+1:end,:) = Aplot;
        elseif direction == "xmax"
            Aplot2 = zeros(size(Aplot,1)*1,size(Aplot,2)*2,size(Aplot,3));
            Aplot2(:,1:size(Aplot,2),:) = Aplot;
            Aplot2(:,size(Aplot,2)+1:end,:) = flip(Aplot,2);
        end
    end

    %% 细化功能
    function Aplot2 = RefineFilter(varargin)
        xPlot=varargin;
        Aplot2 = xPlot{1};
        factor = xPlot{2};

        [nely,nelx,nelz] = size(Aplot2);
        Nelx = nelx*factor;       Nely = nely*factor;     Nelz = nelz*factor;
        %% Elemental information mapping
%         Aplot2(Aplot2 < 0.5) = 0;
%         Aplot2(Aplot2 >= 0.5) = 1;
        Aplot2 = double(reshape(Aplot2,nely,nelx,nelz));
        % xPlot(:,:,1:2) = 1;
        [x,y,z] = ndgrid(1:nely,1:nelx,1:nelz);
        F = griddedInterpolant(x,y,z,Aplot2);
        [X,Y,Z] = ndgrid([1:Nely]./factor,[1:Nelx]./factor,[1:Nelz]./factor);
        xTmpNode0 = F(X,Y,Z);

        rmin = xPlot{3};
        [dy,dx,dz] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
        h = max(0,rmin-sqrt(dx.^2+dy.^2+dz.^2));
        Hs = convn(ones(Nely,Nelx,Nelz),h,'same');
        Aplot2 = convn(xTmpNode0,h,'same')./Hs;
    end
end
end