% =========================================================================
% -- Calculate BER
% -------------------------------------------------------------------------
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [ber] = calculate_ber(par,var,C,s_hat,t,tt)

[~,idxhat] = min(abs(s_hat*ones(1,length(par.symbols))-var.vector*par.symbols).^2,[],2);
bithat = var.bits(idxhat,:);
bits = C.*var.fixed_bits(:,:,t,tt);
ber = sum(bits ~= bithat,2);

end

