function [ massGt ] = mmSLE2Gt( mmSLE )
%MMSLE2GT convert mm of sea level equivalence to ice mass

areaOcean = 362.69*10^6*10^6; %in meters srq
densityIce = 1000; %kg/m^3

masskg = mmSLE/1000*densityIce*areaOcean;
massGt = masskg./10^(12);
end
