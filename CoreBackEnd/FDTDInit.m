function [field,conf,T ] = FDTDInit( conf,delta )
if ~exist('delta','var')
     conf.delta = 3e8./conf.fmax/10;  
else
    conf.delta = delta;
end
conf.deltat = conf.delta/3e8/sqrt(2);    
T = (0:conf.nrOfFrames-1)*conf.deltat;

conf.Xnum = round(conf.x_length/conf.delta);
conf.Ynum = round(conf.y_length/conf.delta);

field.Ez = zeros(conf.Xnum,conf.Ynum);
field.Hx = zeros(conf.Xnum,conf.Ynum-1);
field.Hy = zeros(conf.Xnum-1,conf.Ynum);

field.MuRel   = ones(conf.Xnum*2-1,conf.Ynum*2-1);
field.Sig     = zeros(conf.Xnum*2-1,conf.Ynum*2-1);
field.SigM    = zeros(conf.Xnum*2-1,conf.Ynum*2-1);
field.EpsRel  = ones(conf.Xnum*2-1,conf.Ynum*2-1);

field = repmat(field,1,conf.nrOfFrames);
end

