function [PSD,X,window] = prdgram(sampled_signal,sample_rate,windowtype)
%PRDGRAM Takes a sampled signal and returns the Power Spectral Density
%   estimate using the default parameters of PERIODGRAM, along with a 
%   length normalized X frequency variable, and the window used.
%   
%   Window type can be a default window name such as RECTWIN BARTLETT or
%   BLACKMAN, or a vector of the same length as sampled_signal. Default is
%   RECTWIN
%
%   See also: PERIODOGRAM

defval('windowtype','rectwin')
defval('sample_rate',1)

if ischar(windowtype)
    window = eval(strcat(windowtype,'(length(sampled_signal))'));
else
    window = windowtype;
end
    
[PSD,X] = periodogram(sampled_signal,window,...
    max(256,2^nextpow2(length(sampled_signal))),sample_rate);
end