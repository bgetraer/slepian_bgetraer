function [ anomData, avgMonth ] = anomMonth( data, time)
%ANOMMONTH Normalizes temporally varying data by subtracting the monthly 
% average.
%
% INPUT
%   data    2D or 3D matrix
%   time    date vector
% OUTPUT
%   anomData    2D or 3D matrix of anomalies
%   avgMonth    the monthly mean which has taken out
%
% Last modified: bgetraer@princeton.edu 3/5/2019

data = squeeze(data);

dimen = length(size(data));

switch dimen
    case 3
        if size(data,3)~=length(time)
            error('data and time do not correspond')
        end
        
        if nargin == 3
            avgMonth = zeros(size(data,1),size(data,2),12);
            anomData = zeros(size(data));
            % Take out the average
            for m = 1:12
                forwardD = diff(data,1,3); 
                time = time(1:end-1);
                thisMonth = forwardD(:,:,month(time)==m);
                avgMonth(:,:,m) = mean(thisMonth,3);
                anomData(:,:,month(time)==m) = ...
                    (thisMonth - avgMonth(:,:,m));
            end
        else
            avgMonth = zeros(size(data,1),size(data,2),12);
            anomData = zeros(size(data));
            % Take out the average
            for m = 1:12
                thisMonth = data(:,:,month(time)==m);
                avgMonth(:,:,m) = mean(thisMonth,3);
                anomData(:,:,month(time)==m) = ...
                    (thisMonth - avgMonth(:,:,m));
            end
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

