function [lobeWidth] = LobeWidth(Window,sigLength,Fs)
%WINDOWWIDTHS Summary of this function goes here
%   Detailed explanation goes here

M = sigLength/Fs;
switch Window
    case 'rectwin'
        lobeWidth = 2/M;
    case {'barlett', 'hanning', 'hamming'}
        lobeWidth = 4/M;
    case 'blackmann'
        lobeWidth = 6/M;
    otherwise
        lobeWidth = nan; 
end
end

