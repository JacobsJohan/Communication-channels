function [index] = meter2index(meter,conf)
%meter2index Takes an array with values in meter and gives indexes
%   Takes array of values in meter and returns the indexes using the delta
%   in the configuration structure (for EpsRel and MuRel multiply result by
%   2)
index = round(meter/conf.delta);

end

