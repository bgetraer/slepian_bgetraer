function [] = setworkspace( directory )
%SETWORKSPACE sets the matlab environment 'IFILES' for use with the
%'SLEPIAN_XXXX' function folders.
%
% INPUT
%   directory       the name of the directory containing the 'SLEPIAN_XXXX' 
%                       suite of functions (string)
%
%   sets 'directory' as the working MATLAB environment 'IFILES'
%   adds the 'SLEPIAN_XXXX' suite of functions located in 'directory' to the 
%       current path
%
% Last modified by bgetraer@princeton.edu, 11/23/2017

% default is set for me
if exist('directory','var')==0 
    directory = '/Users/benjamingetraer/Documents/IndependentWork/SH_Workspace';
end

% set IFILES
setenv('IFILES',directory);
% add the slepian directories to the path
addpath(fullfile(getenv('IFILES'),'slepian_alpha'),...
    fullfile(getenv('IFILES'),'slepian_alpha','REGIONS'),...
    fullfile(getenv('IFILES'),'slepian_bravo'),...
    fullfile(getenv('IFILES'),'slepian_delta'),...
    fullfile(getenv('IFILES'),'slepian_echo'),...
    fullfile(getenv('IFILES'),'slepian_foxtrot'),...
    fullfile(getenv('IFILES')));
end