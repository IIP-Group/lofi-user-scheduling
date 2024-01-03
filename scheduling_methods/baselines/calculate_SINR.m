% =========================================================================
% -- Calculate the SINR per user for CSS method
% -------------------------------------------------------------------------
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [SINR_u_final] = calculate_SINR(par,var,H_part)

[~,Us] = size(H_part);

W = pinv(H_part*H_part'+(var.average_N0/var.Es)*eye(par.B))*H_part; % MMSE Equalizer

SINR_u = zeros(Us,1);

for u=1:Us
    des_sig = (var.Es*abs(W(:,u)'*H_part(:,u))^2); % desired signal
    interf = (var.Es*((abs(W(:,u)'*H_part).^2*ones(Us,1))- (abs(W(:,u)'*H_part(:,u))^2))); % interference
    noise = (var.average_N0*norm(W(:,u))^2); % noise

    SINR_u(u) = des_sig/(interf+ noise);
end

SINR_u_final = SINR_u(Us);

end
