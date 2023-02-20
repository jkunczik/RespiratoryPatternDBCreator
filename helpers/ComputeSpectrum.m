function [ P, f, Y, lobeWidth ] = ComputeSpectrum( y, Fs, varargin )
%COMPUTESPECTRUM Summary of this function goes here
%   Detailed explanation goes here


p = inputParser;
addParameter(p,'Window','hamming');
addParameter(p,'MinFreqRes', 0.005);
addParameter(p, 'Truncate', true);
addParameter(p, 'Detrend', true);
p.parse(varargin{:});

Window = p.Results.Window;
L_min = 2^nextpow2(Fs/p.Results.MinFreqRes); 
Truncate = p.Results.Truncate;
Detrend = p.Results.Detrend;

if ~iscell(y)
    y = {y};
    noCellSig = true;
else
    noCellSig = false;
end

numSig = length(y);
Y = cell(1,numSig);
P = cell(1,numSig);

for j=1:numSig
    %% compute fourier transform
    if length(y{j}) > L_min
        L = 2^nextpow2(length(y{j}));   % Make sure the signal is not 
                                        % truncated. The freq. resolution 
                                        % will be better than required
    else
        L = L_min;          % Achieve freq. resolution by padding the 
                            % signal with zeros
    end
%     L = 10e3;
    win = eval([lower(Window) '(' num2str(length(y{j})) ');']);
    sig = y{j}; % repmat(win, 1, size(y{j},2)).*detrend(y{j});
    if Detrend
       sig = detrend(sig); 
    end
    Y{j} = fft(sig, L);
    
    %% Compute auto spectrum
    U = win'*win;
    Sxx = Y{j}.*conj(Y{j})/U; % Auto spectrum.
    
    %% Truncate or shift result
    if Truncate
        % only take the positive side of the spectrum and double it
        if rem(L,2)
          to_select = 1:(L+1)/2;  % ODD
          to_double = 2:length(to_select);
        else
          to_select = 1:L/2+1;    % EVEN
          to_double = 2:length(to_select) -1;
        end
        
        Sxx = Sxx(to_select,:);
        Sxx(to_double,:) = 2*Sxx(to_double,:);
        
        Y{j} = Y{j}(to_select,:);
        Y{j} = Y{j} / length(y{j});
        Y{j}(to_double,:) = 2*Y{j}(to_double,:);
        
        f = Fs*(0:(L/2))/L;
    else
        % shift signal to center it around 0 Hz
        f = (-L/2:L/2-1)*Fs/L;
        Y{j} = fftshift(Y{j});
        Sxx = fftshift(Sxx);
    end
    
    %% compute PSD
    P{j} = Sxx./Fs;
    
end

if noCellSig
   Y = Y{:};
   P = P{:};
end

lobeWidth = LobeWidth(Window,length(y),Fs);

end

