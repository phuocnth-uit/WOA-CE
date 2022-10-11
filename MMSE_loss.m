

function loss = MMSE_loss(Y, Xp, pilot_loc, Nfft, Nps, Nbs, SNR, t_rms, f_max)

 
H_est = MMSE_ideal(Y,Xp,pilot_loc,Nfft,Nps,SNR, t_rms, f_max);

Y_eq = Y./H_est;

Npilot = ceil(Nfft/Nps);
Data_extracted = zeros(Nfft- Npilot, 1);
ip = 0;
for k=1:Nfft
    if mod(k,Nps)==1
        ip=ip+1;
    else
        Data_extracted(k-ip) = Y_eq(k);
    end
end


msg_detected = qamdemod(Data_extracted, 2^Nbs);
Data_ref     = qammod(msg_detected, 2^Nbs);

loss = norm(Data_extracted - Data_ref);

end



