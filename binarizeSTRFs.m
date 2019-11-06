% function [cleanedSTRF]=binarizeSTRFs(inputSTRF)
%
%   FILE NAME   : binarizeSTRFs.m
%   DESCRIPTION : This file converts data points that are > 0.05/# pixels 
%                   as significant (1), else 0
%   
%
% INPUT PARAMS
%   inputSTRF       : original STRF to binarize
%
% RETURNED VARIABLES
%   cleanedSTRF     : binarized version of original STRF
%
% (C) Shannon Lin, Edited Nov 2019

function [cleanedSTRF]=binarizeSTRFs(inputSTRF)
    [rows, cols] = size(inputSTRF);
    numPixels = double(rows * cols);
    threshold = 0.05 / numPixels;
    for i=1:rows
        for j=1:cols
            if(inputSTRF(i,j) > threshold)
                inputSTRF(i,j) = 1;
            else
                inputSTRF(i,j) = 0;
            end
        end
    end
    cleanedSTRF = inputSTRF;
end
