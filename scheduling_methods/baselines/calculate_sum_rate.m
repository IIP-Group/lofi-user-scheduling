% =========================================================================
% -- Calculate the sum-rate for greedy algorithm
% -------------------------------------------------------------------------
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [sum_rate] = calculate_sum_rate(par,var,H_part)

[~,Us] = size(H_part);

W = pinv(H_part*H_part'+(var.average_N0/var.Es)*eye(par.B))*H_part; % MMSE Equalizer

SINR_u = zeros(Us,1);
R_u = zeros(Us,1);

for u=1:Us
    des_sig = (var.Es*abs(W(:,u)'*H_part(:,u))^2); % desired signal
    interf = (var.Es*((abs(W(:,u)'*H_part).^2*ones(Us,1))- (abs(W(:,u)'*H_part(:,u))^2))); % interference
    noise = (var.average_N0*norm(W(:,u))^2); % noise

    SINR_u(u) = des_sig/(interf+ noise);
    R_u(u) = log2(1+real(SINR_u(u))); % rate per user
end

sum_rate = sum(R_u);

end
