% =========================================================================
% -- Fast SURE denoiser
% -------------------------------------------------------------------------
%
% Reference:
%       S. H. Mirfarshbafan, A. Gallyas-Sanhueza, R. Ghods and
%       C. Studer, "Beamspace Channel Estimation for Massive MIMO mmWave
%       Systems: Algorithm and VLSI Design," in IEEE Transactions on Circuits
%       and Systems I: Regular Papers, vol. 67, no. 12, pp. 5482-5495, Dec.
%       2020, doi: 10.1109/TCSI.2020.3023023.
%
% Last Updated: 20/12/2023
%
% -- (c) 2018 Ramina Ghods, Christoph Studer, Seyed Hadi Mirfarshbafan
% -- e-mails: <rghods@cs.cmu.edu, studer@ethz.ch, mirfarshbafan@iis.ee.ethz.ch>
% =========================================================================
function [Hdenoised, mse_avg] = BEACHES(par,Hn, N0, sim_mode)
scalemax = 0.8;
hdenoised = zeros(size(Hn));
SURE = zeros(par.B,1);
mse = zeros(par.U,1);
for uu=1:par.U

    if strcmp(sim_mode, 'antenna_domain')
        hnoisy = fft(Hn(:,uu))/sqrt(par.B);
    else
        hnoisy = Hn(:,uu);
    end

    tau_max = scalemax*max(abs(hnoisy));
    N = length(hnoisy);
    sorth = sort(abs(hnoisy),'ascend');
    recip = 1./sorth;
    recip(sorth < 2^(-10)) = 0;
    cumsum2 = 0;
    cumsuminv = sum(recip);
    tau_opt = inf;
    suremin = inf;
    for bb = 1:N
        tau = sorth(bb);
        SURE(bb) = cumsum2 + (N-bb+1)*tau^2 + N*N0 - 2*N0*(bb-1)-tau*N0*cumsuminv;
        cumsum2 = cumsum2 + sorth(bb).^2;
        cumsuminv = cumsuminv - recip(bb);
        if SURE(bb)<suremin
            suremin = SURE(bb);
            tau_opt = tau;
        end
    end
    tau_opt = min(tau_opt, tau_max);
    hdenoised(:,uu) = exp(1i*angle(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    mse(uu) = suremin/N;
end
if strcmp(sim_mode, 'antenna_domain')
    Hdenoised = ifft(hdenoised)*sqrt(par.B);
else
    Hdenoised = hdenoised;
end
mse_avg = mean(mse);
end