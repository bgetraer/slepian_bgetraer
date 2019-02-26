function [ I, L ] = wavelevelINDEX( C,S )
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
% OUTPUT
%   I       Cell array of coefficient location arrays for each level
%   L       Level array corresponding to the cells in I
%
% SEE ALSO:
%   WAVEDEC2
%
% Last modified: bgetraer@princeton.edu 2/23/2019


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

% Did we lost any coefficients?
if additup~=length(C)
    error("Function error for this decomposition input: lost some coefficients!")
end

end

