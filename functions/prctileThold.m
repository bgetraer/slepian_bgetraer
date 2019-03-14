function [ ptile, level, threshpassindex, CT, S, C, T] = prctileThold( originalimage, target, wname, level )
%PRCTILETHOLD Finds the percentile threshold at which wavelet coefficients
%can be discarded while maintaining the provided target of image invariance.
%   Searches for percentile threshold starting at 0 and increasing by 10,
%   1, 0.1, 0.01, 0.001, 0.0001. and returns the best percentile threshold 
%   maintaining equal or better invariance to the target provided.
%
% INPUT
%   originalimage   the image which is being deconstructed
%   target          the value of bias or invar to be maintained (def. 90%)
%   wname           wavelet to be used in the decomposition (def. haar)
%   level           the level of decomposition (default to highest
%                       possible)
%
%   'test'          outputs a figure demonstrating the WTHCOEF2 function
%
% OUTPUT
%   ptile           the percentile threshold
%   level           the decomposition level used
%   threshpassindex the reference index for eliminating wavelet coeffs
%   CT, S, C, T     wavelet data for recreating the original and 
%                      thresholded images.
%
% SEE ALSO:
%   IMBIAS, IMINVAR, WAVEDEC2, WTHCOEF2, WAVEREC2
%
% last edited by bgetraer@princeton.edu 2/21/2019

defval('target',0.9)
defval('wname','haar')
defval('level',wmaxlev(size(originalimage),wname))

switch nargin
    % Practical implementation
    case {2, 3, 4}
        % initialize ptile to 0th percentile
        ptile = 0;
        % perform the wavelet deconstruction
        [C,S] = wavedec2(originalimage,level,wname);
        
        % absolute values of wavelet coefficients
        abC = abs(C);
        N = 1:level;
        
        n = 6;  % increment go from 10^1 to 10^(-n+2) [ie 10^-4]
        for i = flip(1:n)
            % 10.00 1.00 0.10 0.01
            increment = 10^(i-(n-1));
            invar = 1;   % initialize invariance to 1;
            while invar >= target
                % increment percentile threshold
                ptile = ptile + increment;
                if ptile>100, ptile=100;
                    break
                end
                T = prctile(abC,ptile);
                % Hard thresholding of wavelet coefficients at the threshold level
                CT = wthcoef2('t',C,S,N,repmat(T,level),'h');
                % Reconstruction of the image after thresholding
                testimage = waverec2(CT,S,wname);
                %         fprintf('image invariance = %0.3f \n',iminvar(originalimage, testimage))
                %         fprintf('percentile = %0.3f \n',ptile)
                invar = iminvar(originalimage, testimage);
            end
            % decrement percentile threshold
            ptile = ptile - increment;
        end
        
        % Get the right one
        T = prctile(abC,ptile);
        % Hard thresholding of wavelet coefficients at the threshold level
        CT = wthcoef2('t',C,S,N,repmat(T,level),'h');
        % index for the wavelet coefficients which pass the threshold
        threshpassindex = CT~=0;
        
    case 1
        % Test example of threshold implementation on test image
        if originalimage~='test'
            error('enter full inputs, or "test"');
        end
        % check data directory if files exist
        fileloc = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles/WAVELET_SUPPORT/lena_wavedec.mat';
        if exist(fileloc,'file')
            load(fileloc)
            
            % If the files don't exist, make them!
        else
            fprintf('creating test image files')
            originalimage = imread('lena_std.tif'); % the Lena test image
            originalimage = double(originalimage);
            invarT = [0.1:0.2:0.7 0.8 0.9 0.99];
            wavename = 'haar';
            testimage = zeros(size(originalimage,1),size(originalimage,2),...
                size(originalimage,3),length(invarT));
            
            for i=1:length(invarT) % invariance targets
                for j = 1:size(originalimage,3) % rgb bands
                    [ptile, level] = prctileThold(originalimage(:,:,j), invarT(i),wavename,7);
                    
                    [wdiff,sdiff]=wavedec2(originalimage(:,:,j),level,wavename);
                    abC = abs(wdiff);
                    N = 1:level;
                    T = prctile(abC,ptile);
                    DT = wthcoef2('t',wdiff,sdiff,N,repmat(T,size(N)),'h');
                    testimage(:,:,j,i) = waverec2(DT,sdiff,wavename);
                end
            end
            % change this to the correct directory!!
            fileloc = '/Users/benjamingetraer/Documents/IndependentWork/slepian_bgetraer/datafiles/WAVELET_SUPPORT/lena_wavedec';
            save(fileloc,'originalimage','testimage','invarT');
        end
        % plot
        for i=1:length(invarT) % invariance target
            figure(1)
            clf
            subplot(1,3,1)
            imshow(uint8(originalimage))
            title('Original Image')
            subplot(1,3,2)
            imshow(uint8(testimage(:,:,:,i)))
            title(sprintf('Invar Target: %i%%',round(invarT(i)*100)))
            subplot(1,3,3)
            imshow(double(originalimage) - double(testimage(:,:,:,i)))
            title('Residuals')
            pause(.50)
        end
end
end

