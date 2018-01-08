function [ mmSLE ] = Gt2mmSLE( massGt )
%GT2MMSLE convert ice mass to mm of sea level equivalence

masskg = massGt.*10^(12); %in kg
areaOcean = 362.69*10^6*10^6; %in meters srq
densityIce = 1000; %kg/m^3

mmSLE = 1000*masskg/densityIce./areaOcean; %in mm
end

