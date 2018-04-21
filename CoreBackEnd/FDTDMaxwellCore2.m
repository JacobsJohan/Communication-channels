function [ results ] = FDTDMaxwellCore2(oldvalues,conf,Source,CHXH,CHXE,CHYH,CHYE)
%% Calculation of fields.
    disp([num2str(i),' / ',num2str(conf.numberOfFrames)])
    results.Hx(:,:) =    CHXH.*      field(i).Hx(:,:) -...
                            CHXE.*(     field(i).Ez(:,(1:end-1)+1)  -   field(i).Ez(:,1:end-1)  );
    results.Hy(:,:) =    CHYH.*      field(i).Hy(:,:) +...
                            CHYE.*(     field(i).Ez((1:end-1)+1,:)  -   field(i).Ez(1:end-1,:)  );
                    
    results.Ez(2:end-1,2:end-1) =    CEZE(2:end-1,2:end-1).*     field(i).Ez(2:end-1,2:end-1) +...
                                        CEZH(2:end-1,2:end-1).*(    (field(i+1).Hy(2:end,2:end-1)   -   field(i+1).Hy((2:end)-1,2:end-1))-...
                                                                    (field(i+1).Hx(2:end-1,2:end)   -   field(i+1).Hx(2:end-1,(2:end)-1)));

end

