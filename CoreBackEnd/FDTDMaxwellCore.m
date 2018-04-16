function [ field ] = FDTDMaxwellCore(field,conf,Source)

deltat  = conf.deltat;
delta   = conf.delta;
c       = 3e8;              %Speed of light
mu0     = pi*4e-7;          %Permeabillity of free space 
eps0    = 1/mu0/c/c;        %Permitivity of free space 
Z0      = sqrt(mu0/eps0);   %Free space impedance 
Sc      = c*deltat/delta;  	%Courant number

MUXrel  = field(1).MuRel(1:2:end,2:2:end);
MUYrel  = field(1).MuRel(2:2:end,1:2:end);
SIG     = field(1).Sig(1:2:end,1:2:end);
SIGm    = field(1).SigM(1:2:end,1:2:end);
EPSrel  = field(1).EpsRel(1:2:end,1:2:end);

Chxh = (1-SIGm(:,end-1).*deltat/2./MUXrel)./(1+SIGm(:,end-1).*deltat/2./MUXrel);
Chxe = (Sc/Z0)./MUXrel * 1./(1+SIGm(:,end-1).*deltat/2./MUXrel);
Chyh = (1-SIGm(end-1,:).*deltat/2./MUYrel)./(1+SIGm(end-1,:).*deltat/2./MUYrel);
Chye = (Sc/Z0)./MUYrel * 1./(1+SIGm(end-1,:).*deltat/2./MUYrel);
CEZE = (1-deltat./2*SIG./EPSrel) ./ (1+deltat/2*SIG./EPSrel);
CEZH = 1./(1+deltat/2*SIG./EPSrel) * Z0*Sc./EPSrel;

for i=1:conf.nrOfFrames-1
    disp([num2str(i),' / ',num2str(conf.nrOfFrames)])
    field(i+1).Hx(:,:) =    Chxh.*      field(i).Hx(:,:) -...
                            Chxe.*(     field(i).Ez(:,(1:end-1)+1)  -   field(i).Ez(:,1:end-1)  );
    field(i+1).Hy(:,:) =    Chyh.*      field(i).Hy(:,:) +...
                            Chye.*(     field(i).Ez((1:end-1)+1,:)  -   field(i).Ez(1:end-1,:)  );
                    
    field(i+1).Ez(2:end-1,2:end-1) =    CEZE(2:end-1,2:end-1).*     field(i).Ez(2:end-1,2:end-1) +...
                                        CEZH(2:end-1,2:end-1).*(    (field(i+1).Hy(2:end,2:end-1)   -   field(i+1).Hy((2:end)-1,2:end-1))-...
                                                                    (field(i+1).Hx(2:end-1,2:end)   -   field(i+1).Hx(2:end-1,(2:end)-1)));
   
    for s=1:numel(Source)
       sourceXloc = round(Source(s).X_loc(i)/conf.delta)+1;
       sourceYloc = round(Source(s).Y_loc(i)/conf.delta)+1;
       sourceValue = Source(s).value(i);
%        sourceXloc = max(1,min(conf.Xnum,sourceXloc));
%        sourceYloc = max(1,min(conf.Ynum,sourceYloc));
       field(i+1).Ez(sourceXloc,sourceYloc) = sourceValue;
    end
end

for i = 1:numel(field)
    V(:,:,i)        = field(i).(conf.ToPrint);
    EpsRel(:,:,i)   = field(i).EpsRel;
    MuRel(:,:,i)    = field(i).MuRel;
end

% [M,N,~] = size(V);
% [X,Y]         = meshgrid(     linspace(0,conf.x_length,M),...
%                                 linspace(0,conf.y_length,N));
% [Xq,Yq]      = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
%                                 linspace(0,conf.y_length,conf.Resolution_Y));                
% ToPrintq        = interp2(X,Y,V,Xq,Yq);
% 
% [M,N,~] = size(EpsRel);
% [X2,Y2]         = meshgrid(      linspace(0,conf.x_length,M),...
%                                     linspace(0,conf.y_length,N));
% [Xq2,Yq2]       = meshgrid(     linspace(0,conf.x_length,conf.Resolution_X),...
%                                     linspace(0,conf.y_length,conf.Resolution_Y));
%                     
% EPSrelalpha     = (interp2(X2,Y2,Z2,EpsRel,Xq2,Yq2,Zq2) -1) /   (max(EpsRel(:))-1)/2.5;
% MUrelalpha      = (interp2(X2,Y2,Z2,MuRel,Xq2,Yq2,Zq2)  -1) /   (max(MuRel(:))-1)/2.5;
% 
% temp = ToPrintq(:,:,20:end);
% maxToPrint = max(temp(:));
% minToPrint = min(temp(:));
% absMaxToPrint = max(ToPrintq(:));
% 
% 
% %% Print
% v = VideoWriter(filename,'MPEG-4');
% v.Quality = 100; 
% open(v);
% figure()
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [0, 0, pos(3)*2, pos(4)*2])
% for i = 1:conf.Resolution_T
%     disp(['Frame: ',num2str(i),' / ',num2str(conf.Resolution_T)])
%     surf(Xq(:,:,i),Yq(:,:,i),ToPrintq(:,:,i),...
%             'LineStyle','none',...
%             'FaceColor','interp');
%     hold on 
%     surf(   Xq2(:,:,i),...
%             Yq2(:,:,i),...
%             ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint,...
%             'FaceAlpha','interp',...
%             'AlphaDataMapping','none',...
%             'AlphaData',EPSrelalpha(:,:,1),...
%             'LineStyle','none',...
%             'FaceColor','red');
%     surf(   Xq2(:,:,i),...
%             Yq2(:,:,i),...
%             ones(conf.Resolution_X,conf.Resolution_Y)*absMaxToPrint+0.1,...
%             'FaceAlpha','interp',...
%             'AlphaDataMapping','none',...
%             'AlphaData',MUrelalpha(:,:,1),...
%             'LineStyle','none',...
%             'FaceColor','blue');
%     text(X2(end,end)/10,Y2(end,end)/10,absMaxToPrint+0.2,['time = ',num2str(Zq(1,1,i)),'s']);
%     hold off
%     colorbar;
%     zlim([minToPrint,absMaxToPrint+0.2]);
%     caxis([plotMin,plotMax])
%     view(2)
%     frame = getframe;
%     %img = print('-RGBImage');
%     %frame = im2frame(img);
%     writeVideo(v,frame);
% end
% close(gcf)
% close(v)

end

