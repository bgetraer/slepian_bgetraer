function [ anomData, avgMonth ] = anomMonth( data, time)
%ANOMMONTH Normalizes temporally varying data by subtracting the monthly 
% average.
%
%




dimen = length(size(data));

switch dimen
    case 3
        if size(data,3)~=length(time)
            error('data and time do not correspond')
        end
        avgMonth = zeros(size(data,1),size(data,2),12);
        anomData = zeros(size(data));
        % Take out the average
        for m = 1:12
            thisMonth = data(:,:,month(time)==m);
            avgMonth(:,:,m) = mean(thisMonth,3);
            anomData(:,:,month(time)==m) = ...
                (thisMonth - avgMonth(:,:,m));
        end      
    case 2
        if length(data)~=length(time)
            error('data and time do not correspond')
        end
        if min(size(data))~=1
            error('data must be 1D or 3D')
        end
        avgMonth = zeros(1,12);
        anomData = zeros(size(data));
        % Take out the average
        for m = 1:12
            thisMonth = data(month(time)==m);
            avgMonth(m) = mean(thisMonth);
            anomData(month(time)==m) = ...
                (thisMonth - avgMonth(m));
        end     
end
end

