function [ bool ] = lengthyProcessFlag(description)
%LENGTHYPROCESSFLAG Flag a lengthy process that the user is about to ask
%for and request confirmation. Returns a boolean "true" (continue) or
%"false" (skip)
%
% USEAGE:
%   if LENGTHYPROCESSFLAG
%       expensiveFunction();
%   end
%
% Last modified by bgetraer@princeton.edu 3/12/2019

switch nargin
    case 0
        prompt = 'Lengthy Process has been flagged. Continue? (y/n)\n';
    case   1
        prompt = sprintf('Lengthy Process "%s" has been flagged. Continue? (y/n)\n',...
            description);
end

clc
command = input(prompt,'s');

if strcmp(command,'y') || strcmp(command,'yes') || ...
        strcmp(command,'Y') || strcmp(command,'Yes')
    bool = 1;
elseif strcmp(command,'n') || strcmp(command,'no') || ...
        strcmp(command,'N') || strcmp(command,'No')
    bool = 0;
else
    error('Invalid input!')
end
end

