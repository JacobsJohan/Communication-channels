function [field,conf,T ] = FDTDInit_fastfading( conf,delta )
if ~exist('delta','var') %Calculate spatial delta if not given in function.
     conf.delta = 3e8./conf.fmax/10;  
else
    conf.delta = delta;
end
conf.deltat = conf.delta/3e8/sqrt(2);   %Calculate temporal delta
T = (0:conf.nrOfFrames-1)*conf.deltat;  %Create temporal axis

conf.x_length       = 500*conf.delta;
conf.y_length       = 500*conf.delta;

conf.Xnum = round(conf.x_length/conf.delta);
conf.Ynum = round(conf.y_length/conf.delta);

%Initializing fields
field.Ez = zeros(conf.Xnum,conf.Ynum);
field.Hx = zeros(conf.Xnum,conf.Ynum-1);
field.Hy = zeros(conf.Xnum-1,conf.Ynum);

%Initialize properties of space
field.MuRel   = ones(conf.Xnum*2-1,conf.Ynum*2-1);
field.Sig     = zeros(conf.Xnum*2-1,conf.Ynum*2-1);
field.SigM    = zeros(conf.Xnum*2-1,conf.Ynum*2-1);
field.EpsRel  = ones(conf.Xnum*2-1,conf.Ynum*2-1);

end

