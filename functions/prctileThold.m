function [ ptile, level ] = prctileThold( originalimage, target, wname, level )
%PRCTILETHOLD Finds the percentile threshold at which wavelet coefficients
%can be discarded while maintaining the provided target of image invariance.
%   Searches for percentile threshold starting at 0 and increasing by 10,
%   1, 0.1, and 0.01 and returns the best percentile threshold maintaining
%   equal or better invariance to the target provided.
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
%
% SEE ALSO:
%   IMBIAS, IMINVAR, WAVEDEC2, WTHCOEF2, WAVEREC2
%
% last edited by bgetraer@princeton.edu 2/17/2019

defval('target',0.9)
defval('wname','haar')
defval('level',floor(log(min(size(originalimage)))/log(2)))

switch nargin
    % Practical implementation
    case {2, 3, 4}
        % initialize ptile to 0th percentile
        ptile = 0;
        
        % perform the wavelet deconstruction
        [C,S] = wavedec2(originalimage,level,wname);
        
        % absolute values of wavelet coefficients
        abwdiff = abs(C);
        N = 1:level;
        
        % initialize test image
        testimage = originalimage;
        
        % loop to find percentile threshold
        for i = flip(1:4)
            % 10.00 1.00 0.10 0.01
            increment = 10^(i-3);
            
            while iminvar(originalimage, testimage) >= target
                % increment percentile threshold
                ptile = ptile + increment;
                T = prctile(abwdiff,ptile);
                % Hard thresholding of wavelet coefficients at the threshold level
                NC = wthcoef2('t',C,S,N,repmat(T,level),'h');
                % Reconstruction of the image after thresholding
                testimage = waverec2(NC,S,wname);
                %         fprintf('image invariance = %0.3f \n',iminvar(originalimage, testimage))
                %         fprintf('percentile = %0.3f \n',ptile)
            end
            % decrement percentile threshold
            ptile = ptile - increment;
            % reset test image
            testimage = originalimage;
        end
        
    % Test example of threshold implementation
    case 1
        if originalimage~='test'
            error('enter full inputs, or "test"');
        end
        % load image
        originalimage = imread('lena_std.tif');
        originalimage_recon = originalimage;
        wname = 'haar';
        ptile = 99;
        level = floor(log(min(size(originalimage(:,:,1)))/log(2)));
        
        % threshold each rgb band
        for i = 1:3
            % perform the wavelet deconstruction
            [C,S] = wavedec2(originalimage(:,:,i),level,wname);
            % absolute values of wavelet coefficients
            abwdiff = abs(C);
            N = 1:level;
            T = prctile(abwdiff,ptile);
            % Hard thresholding of wavelet coefficients at the threshold level
            NC = wthcoef2('t',C,S,N,repmat(T,level),'h');
            % Reconstruction of the image after thresholding
            originalimage_recon(:,:,i) = waverec2(NC,S,wname);
        end
        
        % plot
        figure(1)
        clf
        subplot(1,3,1)
        imshow(originalimage)
        title('Original Image')
        colormap('bone')
        subplot(1,3,2)
        imshow(originalimage_recon)
        title(sprintf('WaveRec Thresholded at %ith percentile',...
            ptile))
        subplot(1,3,3)
        imshow(double(originalimage) - double(originalimage_recon))
        title('Residuals')
end
end

