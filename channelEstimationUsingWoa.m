

% Channel estimation using LS, WOA and MMSE 
% Number of OFDM Symbol = 1e1
% Channel model: TDLC-300

clc, clear; close all;
methods = {'LS ', 'WOA', 'MMSE'}            % Channel estimation methods

snrRange            = 0:5:30;               % Signal to noise ratio in dB
numSymbol           = 1e1;                  % Number of symbols
numFft              = 4096;                 % Size of DFT 
numCp               = numFft/4;             % Number of CP
subCarrierSpacing   = 30e3;                 % Subcarrier Spacing
numBitPerSym        = 4;                    % Number of bits per (modulated) symbol
numSymPerPilot      = 12;                   % Number of (modulated) symbol per pilots
numBitBerSecond     = 1e3;                  % Number of bits per second
signalEnergy        = 100;                  % Energy of signal

% Propagation Channel Model Configuration
% Create a TDL channel model object and specify its propagation characteristics.
numTapEst   = 400;                          % Number of est. channel taps
numTap      = 320;                          % Number of true channel taps

% TDLC300-100
tapDelay    = [0    65 70   190  195  200  245  325  520  1045  1510  2595];    % in ns
tapPower    = [-6.9 0  -7.7 -2.5 -2.4 -9.9 -8.0 -6.6 -7.1 -13.0 -14.2 -16.0];   % in dB 

% WOA Alg
maxIter     = 8;                            % maximum number of generations
numAgent    = 8;                            % Number of search agents

ub          = [50   100 400];               % [ub1,ub2,...,ubn] where ubn is the upper bound of variable n
lb          = [0    20  0];                 % [lb1,lb2,...,lbn] where lbn is the lower bound of variable n
dim         = 3;                            % Number of variables
positions   = rand(numAgent, dim).*(ub-lb) + lb;

visualization = 0;
saveOrNot = 0;                              % = 1 for save

sampRate    = numFft*subCarrierSpacing;     % Sample rate
numPilot    = ceil(numFft/numSymPerPilot);  % Number of pilots per OFDM symbol
pilotLoc    = zeros(numPilot, 1);           % Pilot's Location

pathLoss    = zeros(numTap, 1);         
tapSample   = round(tapDelay*1e-9*sampRate);  
pathLoss(tapSample+1) = 10.^(tapPower/10);  % Path loss of channel

M           = 2^numBitPerSym;               % M - QAM
A           = sqrt(3/2/(M-1)*signalEnergy); % QAM normalization factor
Nofdm       = numFft + numCp;               % Number of OFDM
numData     = numFft - numPilot;            % Number of data

MSEs_snr    = zeros(length(snrRange),length(methods));
ber_snr     = zeros(length(snrRange),length(methods));
fileIdx     = getFileId(saveOrNot);
tic 

for snrIdx = 1:length(snrRange)
    SNR = snrRange(snrIdx);
    er = zeros(1,length(methods));
    MSE = zeros(1,length(methods));
    for nsym=1:numSymbol

        msgint = randi([0 M-1],numFft-numPilot,1);      % Symbol generation
        data = qammod(msgint, M);


        % Add pilot
        p = randi([0, M-1], numPilot, 1);               % Pilot sequence generation
        pilot = qammod(p, M);
        ip = 0;
        X = zeros(numFft, 1);

        for k=1:numFft
            if rem(k,numSymPerPilot)== 1
                ip = ip+1;
                X(k)=pilot(floor(k/numSymPerPilot)+1);  % For pilot
                pilotLoc(ip) = k;                       % For pilot location
            else
                X(k) = data(k-ip);                      % For data
            end
        end

        % OFDM
        x = ifft(X,numFft);                             % IFFT
        xt = [x(numFft-numCp+1:numFft); x];             % add CP

        % PA
        tx = A*xt;
        signalPowerdB = 10*log10(cov(tx));

        % Channal gain
        h = (randn(numTap, 1)+1j*randn(numTap, 1))...
            .*sqrt(pathLoss/2);                        % Channel gain
        H = fft(h,numFft);                             % True channel frequency respond
        H_power_dB = 10*log10(abs(H.*conj(H)));        % True channel power in dB
        y_channel = conv(tx,h);                        % Channel path (convolution)

        % Add noise
        rx = awgn(y_channel, SNR, 'measured');
        % rx = y_channel + 1/(sqrt(2.0)*10^(SNR/20))*complex(randn(size(y_channel)),randn(size(y_channel)));

        % sto = sto_est(rx, numFft, numCp);
        % Receiver
        y = rx(numCp+1:Nofdm);                         % Remove CP
        Y = fft(y);                                    % FFT


        % Channel estimation
        for methodIdx = 1:length(methods)
            method = methods{methodIdx};
            if     method(1) == 'L'
                % LS estimation with linear interpolation
                H_est = LS_CE(Y,pilot,pilotLoc,numFft, 'linear');

            elseif method(1) == 'W'
                % WOA estimation
                [H_est, positions] = WOA_CE(Y,pilot,pilotLoc,numFft, numSymPerPilot, numBitPerSym, ...
                        positions,  numAgent, maxIter, lb, ub, dim);
          
            elseif method(1) == 'M'
                % MMSE estimation
                H_est = MMSE_CE(Y,pilot,pilotLoc,numFft,numSymPerPilot,h,SNR);
            end


            if method(end) == 'T'
                h_est = ifft(H_est);                   % Esti channel gain
                h_est = h_est(1:numTapEst);            % N-tap channel gain
                H_est = fft(h_est,numFft);             % DFT-based channel estimation
            end

            H_est_power_dB = ...
                10*log10(abs(H_est.*conj(H_est)));     % Esti channel power in dB


            Y_eq = Y./H_est;
            Data_extracted = zeros(numFft- numPilot, 1);
            ip = 0;

            for k=1:numFft
                if mod(k,numSymPerPilot)==1
                    ip=ip+1;
                else
                    Data_extracted(k-ip)=Y_eq(k);
                end
            end

            msg_detected = qamdemod(Data_extracted, M);
            bitDetected = de2bi(msg_detected, numBitPerSym);
            bitTrans    = de2bi(msgint, numBitPerSym);

            er(methodIdx) = er(methodIdx) + sum(sum(bitDetected~=bitTrans));

            MSE(methodIdx) = MSE(methodIdx) + (H-H_est/A)'*(H-H_est/A);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if visualization
                figure(snrIdx)
                subplot(1, length(methods), methodIdx)
                hold on
                scatter(real(Data_extracted), imag(Data_extracted), 'b.')
                title(methods{methodIdx})
                pause(0)
            end

        end
    end

    MSEs = MSE/(numFft*numSymbol);
    MSEs_snr(snrIdx, :) = MSEs;
    ber_snr(snrIdx, :) = er/(numSymbol*numData*numBitPerSym);
    
  
    for methodIdx = 1:length(methods)
        str = sprintf('SNR = %2.0f dB: BER of %s \t= %6.5f\n', SNR, methods{methodIdx}, ber_snr(snrIdx, methodIdx));
        fprintf(str)
        if saveOrNot
            if methodIdx == 1

            fprintf(fileIdx, '\n-----------------------------------------\n');
            end

            fprintf(fileIdx, str);
        end
    end
end

if saveOrNot  
    fclose(fileIdx);
end

figure
subplot(121)
semilogy(snrRange, ber_snr)
legend(methods{:})
xlabel('SNR')
ylabel('BER')
grid on

subplot(122)
semilogy(snrRange, MSEs_snr)
legend(methods{:})
xlabel('SNR')
ylabel('MSE')
grid on
toc

% Local Functions
function [H_interpolated] = interpolate(H,pilot_loc,Nfft,method)
% Input: H = Channel estimate using pilot sequence
% pilot_loc = Location of pilot sequence
% Nfft = FFT size
% method = ’linear’/’spline’
% Output: H_interpolated = interpolated channel
if pilot_loc(1)>1
    slope = (H(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
    H = [H(1)-slope*(pilot_loc(1)-1); H]; pilot_loc = [1; pilot_loc];
end

if pilot_loc(end) <Nfft
    slope = (H(end)-H(end-1))/(pilot_loc(end)-pilot_loc(end-1));
    H = [H; H(end)+slope*(Nfft-pilot_loc(end))];
    pilot_loc = [pilot_loc; Nfft];
end

if lower(method(1))=='l'
    H_interpolated = interp1(pilot_loc,H,[1:Nfft]', 'linear');
else
    H_interpolated = interp1(pilot_loc,H,[1:Nfft]', 'spline');
end
end

function H_LS = LS_CE(Y,Xp,pilot_loc,Nfft,int_opt)
% LS channel estimation function
% Inputs:
% Y = Frequency-domain received signal
% Xp = Pilot signal
% pilot_loc = Pilot location
% N = FFT size
% Nps = Pilot spacing
% int_opt = ’linear’ or ’spline’
% output:
% H_LS = LS Channel estimate

LS_est = Y(pilot_loc)./Xp; % LS channel estimation
if lower(int_opt(1))=='l'
    method='linear';
else
    method='spline';
end
% Linear/Spline interpolation
H_LS = interpolate(LS_est,pilot_loc,Nfft,method);
end


function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
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

% Calculate RMS delay spread
Ph = h.*conj(h);
Ptotal = h'*h;
t_sym = 1*(0:length(h)-1)';
t_mean = sum(t_sym.*Ph/Ptotal);
t_cov = sum(t_sym.^2.*Ph/Ptotal);
t_rms = sqrt(t_cov-t_mean^2);
f_max = 100;

H_MMSE = MMSE_ideal(Y,Xp,pilot_loc,Nfft,Nps,SNR, t_rms, f_max);
end


function [H_WOA, Positions] = WOA_CE(Y,Xp,pilot_loc,Nfft,Nps, Nbs, ...
                 Positions, NumAgent, Max_iter, lb, ub, dim)    

% fobj              = @CostFunction
% dim               = number of your variables
% Max_iteration     = maximum number of generations
% SearchAgents_no   = number of search agents
% lb                = [lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub                = [ub1,ub2,...,ubn] where ubn is the upper bound of variable n


fobj = @ (x) MMSE_loss(Y, Xp, pilot_loc, Nfft, Nps, Nbs, x(1), x(2), x(3));

x = WhaleOptAlg(NumAgent,Max_iter,lb,ub,dim,fobj);
SNR   = x(1);
t_rms = x(2);
f_max = x(3);

H_WOA = MMSE_ideal(Y,Xp,pilot_loc,Nfft,Nps,SNR, t_rms, f_max);

end

function fileId = getFileId(enable)
if enable == 0
    fileId = 0;
    return 
end

ctime       = clock;
cmonth     = ctime(2); smonth  = num2str(cmonth);
cday       = ctime(3); sday    = num2str(cday);
chour      = ctime(4); shour   = num2str(chour);
cminute    = ctime(5); sminute = num2str(cminute);

fileName = ['CE', smonth, sday, shour, sminute, '.txt'];

fileId = fopen(fileName, 'w');
end



