%Copyright (c) 2011, Universidad Industrial de Santander, Colombia
%University of Delaware
%All rights reserved.
%@author: Sergio Pino
%@author: Henry Arguello
%Website: http://www.eecis.udel.edu/
%emails  : sergiop@udel.edu - henarfu@udel.edu
%Date   : Sept, 2011

function [ fading ] = fadingEstimator( sine )
% FADINGESTIMATOR estimate fading, asuming fading cyclostationary
%   Estimate alpha from a sine wave
% input:
%   - sine = sine signal
% output:
%   - fading = estimated value of the fading

    %% ALPHA estimate alpha, asume alpha cyclostationary

    % Estimate alpha from a sine wave
    sineFix = abs(sine(ceil(1/3*end):ceil(2/3*end)));
    % floor((100e6/512)/(2*1000)) since we generate a sine with freq 1000Hz
    % and it is sampled at 100e6/512 Hz
    fading = mean(findpeaks(sineFix, 'minpeakdistance', floor((100e6/512)/(2*1000))));

end

