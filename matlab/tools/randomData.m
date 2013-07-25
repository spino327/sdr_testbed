%% 2011 - Sergio Pino
%   University of Delaware

function [a_k, mx, var, p0, p1] = randomData(len)
%randomData Generate random bits for test the BPSK implementation
%   in: - len = number of desired bits
%   out: - random bits

%     a_k = randint(len,1);
    a_k = randi([0, 1], 1, len);

    mx = mean(a_k);
    vx = var(a_k);
    fprintf('mx = %f, var = %f \n', mx, vx);

    p0 = sum(a_k == 0)/len;
    p1 = sum(a_k == 1)/len;

    fprintf('p(0) = %d \n', p0);
    fprintf('p(1) = %d \n', p1);

end