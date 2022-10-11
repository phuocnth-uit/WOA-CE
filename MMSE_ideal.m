

function H_MMSE = MMSE_ideal(Y,Xp,pilot_loc,Nfft,Nps,SNR, t_rms, f_max)
% MMSE channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% Nfft = FFT size
% Nps = Pilot spacing
% h = Channel impulse response
% SNR = Signal-to-Noise Ratio[dB]
% output:
% H_MMSE = MMSE channel estimate

Nsym    = 1;
Npilot  = ceil(Nfft/Nps);
df      = 1/Nfft;
kr      = -Nfft+1 : Nfft-1;
lr      = -Nsym+1 : Nsym-1;

snr = 10^(SNR/10);
H_LS = Y(pilot_loc,1)./Xp;      % LS estimate


% rf[k - k'] and  rt[l - l']
rf = 1./(1 + 1j*2*pi*t_rms*df*kr);
rt = besselj(0, 2*pi*f_max*Nfft*lr);

% Cross-corr matrix
Rpp = zeros(Npilot);
for i = 1:Npilot
    for j = 1:Npilot
        k = pilot_loc(i)-pilot_loc(j) + Nfft;
        l = 0 + Nsym;

        Rpp(i,j) = rf(k)*rt(l);
    end
end

Rhp = zeros(Nfft, Npilot);
for i = 1:Nfft
    for j = 1:Npilot
        k = i-pilot_loc(j) + Nfft;
        l = 0 + Nsym;

        Rhp(i,j) = rf(k)*rt(l);
    end
end

H_MMSE = Rhp/(Rpp+eye(Npilot)/snr)*H_LS;

end
