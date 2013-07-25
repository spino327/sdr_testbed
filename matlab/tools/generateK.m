function [ k ] = generateK( csnr, performance)
%GENERATEK Calculate the factor k, given the csnr and the performance of
%the system for the model y = k*alpha*st + nt
%   [ k ] = generateK( csnr, performance)
%   in: - csnr = target chanel SNR in dB
%       - performance = system performance when k = 1 not in dB, current
%       csnr

    k = sqrt( exp(csnr*log(10)/10) / performance );

end

