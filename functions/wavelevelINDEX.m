function [ I, L, coefinlevel ] = wavelevelINDEX( C,S )
%WAVELEVELINDEX Returns the location of the coefficients of the various
% decomposition levels within the coefficient array C.
%       Coefficients are stored in array C in order by increasing level of
%   spatial resolution. For a Haar wavelet decomposition of max level, the
%   first coefficient
%       The resolution of each level is stored in S, starting with the
%   dimension of the lowest resolution "approximation coefficients" (which
%   have a single value over their support), followed by the dimension of
%   the  "detail coefficients" (which vary horizonatally, vertically, and
%   diagonally), and followed finally by the dimension of the original
%   image.
%       Thus in the C array, there exist one set of approximation
%   coefficients for the lowest resolution level, and 3 sets of of
%   detail coefficients for each resolution level.
%
% INPUT
%   C       The coefficient array from wavedec2
%   S       The size array from wavedec2
%
%   'test'          outputs a figure demonstrating the WTHCOEF2 function
%
% OUTPUT
%   I       Cell array of coefficient location arrays for each level
%   L       Level array corresponding to the cells in I
%
% SEE ALSO:
%   WAVEDEC2
%
% Last modified: bgetraer@princeton.edu 2/23/2019

switch nargin
    case 2
        % extract the number of coefficients in each level from the size array
        coefinlevel = flip(prod(S(2:end-1,:),2));
        % the levels in order of high to low resolution
        L = 1:1:(length(coefinlevel));
        % the index of level coefficient locations within C
        I = cell(1,length(coefinlevel));
        
        additup = 0;    % track the total number of coefficients
        % loop over each level
        for i = flip(L)
            if i == length(coefinlevel)
                % lowest resolution level has approximation and detail coefficients
                I{i} = 1: coefinlevel(i) + 3*coefinlevel(i);
            else
                % all other resolution levels only have detail coefficients
                I{i} = I{i+1}(end) + (1:3*coefinlevel(i));
            end
            additup = additup + length(I{i});
        end
        
        % Did we lost any coefficients? (maybe for other wavelet types)
        if additup~=length(C)
            error("Function error for this decomposition input: lost some coefficients!")
        end
        
    case 1
        % Test example of threshold implementation on test image
        if C~='test'
            error('enter full inputs, or "test"');
        end
        
        clear C
        
        originalimage = imread('lena_std.tif'); % the Lena test image
        originalimage = double(originalimage(:,:,3));
        
        sz = size(originalimage);
        wavename = 'haar';
        level = wmaxlev(sz,wavename);
        [C,S]=wavedec2(originalimage,level,wavename);
        Cbottomup = C*0;
        Ctopdown = C*0;
        
        [ I, L ] = wavelevelINDEX( C,S );
        
        clf
        for i=1:length(L)
            clf
            hold on
            subplot(1,2,1)
            Cbottomup(I{L(i)}) = C(I{L(i)});
            imj = waverec2(Cbottomup,S,wavename);
            imshow(imj,[])
            L(i);
            ttl = "Levels 1:%i";
            title(sprintf(ttl,L(i)))
            text(.5*size(imj,1),1.05*size(imj,2),'Bottom-up Reconstruction',...
                'horizontalalignment','center')
            text(1.15*size(imj,1),-.25*size(imj,2),'\bf Haar Wavelet Decomposition',...
                'horizontalalignment','center','fontsize',14)
            
            subplot(1,2,2)
            j = length(L)+1-i;
            Ctopdown(I{L(j)}) = C(I{L(j)});
            imj = waverec2(Ctopdown,S,wavename);
            
            imshow(imj,[])
            L(j);
            ttl = "Levels %i:%i";
            title(sprintf(ttl,L(j),length(L)))
            text(.5*size(imj,1),1.05*size(imj,2),'Top-down Reconstruction',...
                'horizontalalignment','center')
            pause(1)
        end
        
end

