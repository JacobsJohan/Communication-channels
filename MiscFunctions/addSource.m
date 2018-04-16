function [ Sources ] = addSource( Sources,conf,X_loc,Y_loc,frequency,value )
if frequency>conf.fmax
    error('Max Frequency too small')
end
if numel(X_loc) == 1
    X_loc = X_loc * ones(1,numel(value)); % Makes a vector with all x indices (in meter) of the source for each timeframe.
end
if numel(Y_loc) == 1
    Y_loc = Y_loc * ones(1,numel(value));   % Makes a vector with all y indices (in meter) of the source for each timeframe.
end
if ~isfield(Sources,'X_loc') % isfield True if field is in structure array.
    Sources(1).X_loc = X_loc;
    Sources(1).Y_loc = Y_loc;
    Sources(1).value = value;
else
    Sources(end+1).X_loc = X_loc;
    Sources(end).Y_loc = Y_loc;
    Sources(end).value = value;
end
end

