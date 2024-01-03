% =========================================================================
% -- Generates channel with power control
% -------------------------------------------------------------------------
% Last Updated: 20/12/2023
%
% -- (c) 2023 Victoria Palhares, Christoph Studer
% -- e-mails: <palhares@iis.ee.ethz.ch, studer@ethz.ch>
% =========================================================================

function [var] = gen_channel(par,var)

%% Obtaining channel
FileLoad='Channels_B_16_U_16_Ch_100.mat'; % Channels obtained on Wireless Insite
load (FileLoad,'H');

var.G = H(:,1:par.U,var.ch);

% Power control
if (strcmp(par.pwr_cntrl,'on')) % with power control
    gpwr= sum(abs(var.G).^2,1); % square of L2 norm of each column of G
    gmin = min(gpwr); % minimum of all norms
    p = 10^(par.P_db/10);
    lambda = min([gpwr;gmin*p*ones(size(gpwr))])./gpwr;
    Lambda_vector = sqrt(lambda);
elseif (strcmp(par.pwr_cntrl,'off')) % without power control
    Lambda_vector = ones(par.U,1);
end

var.H = var.G*diag(Lambda_vector); % Applying power control matrix

% Checking range between weakest and strongest user
gpwr= sum(abs(var.H).^2,1);
gpwr_db = 10.*log10(gpwr); % [dB]
dynamic_range = max(gpwr_db)-min(gpwr_db); % [dB]

end