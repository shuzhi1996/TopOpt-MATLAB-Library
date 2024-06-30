classdef OpenMP_Plot2D
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
% class for Model
% How to use:
%   Type SetModel.[function name] on MATLAB command window. 
%-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

methods (Static)
    %/////////////////////////////// class function ////////////////////////////////
    function PrintingPlot(xPhysi,phie,NumPrint,type,FlagBar,FlagTitle,varargin)
        [Nely,Nelx] = size(xPhysi);
        phimax = max(phie(:));  phie(xPhysi<0.45)=NaN;
        [X,Y] = meshgrid(1:Nelx,1:Nely);
        [C,h] = contourf(X,Y,flipud(phie/phimax),NumPrint); axis equal; axis tight; axis off; hold on;
        colormap(type);

        if FlagBar > 0
            colorbar('eastoutside');
        end

        if FlagTitle > 0
            title(strcat('Itr:  ',num2str(varargin{1}),'  Obj:  ',num2str(round(varargin{2},3)),...
                '  Thick:  ',num2str(round(varargin{3},3)),'  Div:  ',num2str(round(varargin{4},3))))
        end
    end

    function PlotNode(TedofMat,xPhysi,U,scale,max)
        [Nely,Nelx] = size(xPhysi);
        [X,Y] = meshgrid(0:Nelx,Nely:-1:0);
   
        U = U*scale;

        clear Sa
        gcf22 = figure(22);   clf
        set(gcf22,'name','Results for The Lattice','numbertitle','off','color','w');
        X1 = X(:) + U(1:2:end-1,1);
        Y1 = Y(:) + U(2:2:end,1);
        Sa.Vertices = [X1(:),Y1(:)];
        Sa.Faces = TedofMat(xPhysi(:)>0.5,:);
        Sa.FaceVertexCData = sqrt(U(1:2:end-1).^2 + U(2:2:end).^2);
        Sa.FaceColor = 'interp';
        Sa.EdgeColor = 'k';
        Sa.LineStyle ='none';
        patch(Sa)
        axis off equal; hold on;
        colormap(jet); colorbar('southoutside');clim([0 max*scale])
    end

    function PlotEle(eta,angle1,angle2,varargin)
        xPlot=varargin;
        Aplot2 = xPlot{1};

        % Aplot2(1:1,1:1,1:1)=1;
        % Aplot2(end,end,end)=1;
        % Bplot = permute(Aplot2,[1 3 2]);

        isovals = Aplot2;%shiftdim(Bplot, 2);
        isovals = smooth3(isovals,'box',1);

        [F1,V1] = isosurface(isovals,eta);
        [F2,V2] = isocaps(isovals,eta);
        F3 = [F1;F2+size(V1,1)];
        V3 = [V1;V2];

        % Sa = isosurface(isovals,eta);
        Sa.Vertices = V3;
        Sa.Faces = F3;

        SS = patch(Sa);
        Bplot22 = xPlot{2};%permute(xPlot{2},[1 3 2]);
        % Bplot22 = shiftdim(Bplot22, 2);
        isocolors(Bplot22,SS);
        SS.FaceColor = 'interp';
        SS.EdgeColor = 'none';

        view([angle1,angle2]); axis equal tight off; camlight; drawnow
        colormap(xPlot{3}); 
%         colorbar('eastoutside');
        ax = gca;
        ax.Interactions = [rotateInteraction dataTipInteraction];
    end
    function Aplot2 = flipFilter(Aplot)
        % xz
        Aplot2 = zeros(size(Aplot,1)*2,size(Aplot,2)*1,size(Aplot,3));
        Aplot2(1:size(Aplot,1),1:size(Aplot,2),:) = flipud(Aplot);
        Aplot2(size(Aplot,1)+1:end,1:size(Aplot,2),:) = Aplot;
        % Aplot2(:,1:size(Aplot,2),:) = flip(Aplot2(:,1:size(Aplot,2),:),2);
        % Aplot2(:,size(Aplot,2)+1:end,:) = flip(Aplot2(:,1:size(Aplot,2),:),2);
    end

    function Plot(loop,ObjectRec,NumCons,varargin)
        ConsRec = varargin;
        %% Plotting Data
        yyaxis left;
        plot(1:loop,ObjectRec);
        ylabel('Objective');
        yyaxis right
        for mm = 1:NumCons
        plot(1:loop,ConsRec{mm}); hold on
        end
        ylabel('Constraint');
        xlabel('Iteration');grid on;
    end
end
end