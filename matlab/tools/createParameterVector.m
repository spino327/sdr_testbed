
function [ vector ] = createParameterVector( position, varargin )
%CREATEPARAMETERVECTOR Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 2
       return 
    end
    
    vector = zeros(1, nargin - 1);
    
    for i = 1:length(vector)
        tmp = varargin(i);
        tmp = tmp{1};
        vector(i) = tmp(position);
    end

end