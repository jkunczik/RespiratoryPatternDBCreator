function [res] = ReferenceSig(fs_disp, plot)
    if nargin < 2
        plot = false;
    end
    
    fs = 30.35;
    t = 0:1/fs:60-1/fs;

    sig = [];
    label = [];
    rr = [];
    a = [];
    
    % 0 initial gap (5 s)
    AddGap(5);

    % 1 bradypnea (5 BPM)
    sig = [sig, sin(2*pi*5/60*t)];
    AddLabel(length(t), Pattern.BRADYPNEA);
    AddFeatures(length(t), 5, 1);
    AddGap(6);

    % 2 eupnea (12 BPM)
    sig = [sig sin(2*pi*12/60*t)];
    AddLabel(length(t), Pattern.EUPNEA);
    AddFeatures(length(t), 12, 1);
    AddGap(5.25);

    % 3 hypopnea (12 BPM)
    sig = [sig 0.25*sin(2*pi*12/60*t)];
    AddLabel(length(t), Pattern.HYPOPNEA);
    AddFeatures(length(t), 12, 0.25);
    AddGap(6);

    % 4 hyperpnea (10 BPM)
    sig = [sig 2.5*sin(2*pi*10/60*t)];
    AddLabel(length(t), Pattern.HYPERPNEA);
    AddFeatures(length(t), 10, 2.5);
    AddGap(5.75);

    % 5 tachypnea (35 BPM)
    sig = [sig sin(2*pi*35/60*t)];
    AddLabel(length(t), Pattern.TACHYPNEA);
    AddFeatures(length(t), 35, 1);
    AddGap(5.25);

    % 6 eupnea (15 BPM)
    sig = [sig sin(2*pi*15/60*t)];
    AddLabel(length(t), Pattern.EUPNEA);
    AddFeatures(length(t), 15, 1);
    AddGap(5.75);

    % 7 kussmaul (30 BPM)
    sig = [sig 2.5*sin(2*pi*30/60*t)];
    AddLabel(length(t), Pattern.KUSSMAUL);
    AddFeatures(length(t), 30, 2.5);
    AddGap(5.5);

    % 8 cheyne stokes (20 BPM)
    t = 0:1/fs:15-1/fs;
    cheyne = 4.5/15*t.*sin(2*pi*20/60*t);
    cheyne_rr = ones(1, length(t))*20;
    cheyne_a = 4.5/15*t;
    cheyne = [cheyne zeros(1, round(20.25*fs))];
    cheyne_rr = [cheyne_rr zeros(1, round(20.25*fs))];
    cheyne_a = [cheyne_a zeros(1, round(20.25*fs))];
    cheyne = [cheyne, cheyne];
    cheyne_rr = [cheyne_rr, cheyne_rr];
    cheyne_a = [cheyne_a, cheyne_a];
    sig = [sig cheyne];
    rr = [rr cheyne_rr];
    a = [a cheyne_a];
    AddLabel(length(cheyne), Pattern.CHEYNE_STOKES);
    AddGap(5.25);
    

    % 9 tachypnea (35 BPM)
    t = 0:1/fs:60-1/fs;
    sig = [sig 0.25*sin(2*pi*35/60*t)];
    AddLabel(length(t), Pattern.TACHYPNEA);
    AddFeatures(length(t), 35, 0.25);
    AddGap(5.25);

    % 10 biot (15 BPM)
    t = 0:1/fs:20-1/fs;
    biot = sin(2*pi*15/60*t);
    biot_rr = ones(1, length(t))*15;
    biot_a = ones(1, length(t));
    biot = [biot zeros(1, round(20.25*fs))];
    biot_rr = [biot_rr zeros(1, round(20.25*fs))];
    biot_a = [biot_a zeros(1, round(20.25*fs))];
    biot = [biot, biot];
    biot_rr = [biot_rr, biot_rr];
    biot_a = [biot_a, biot_a];
    sig = [sig biot];
    rr = [rr biot_rr];
    a = [a biot_a];
    AddLabel(length(biot), Pattern.BIOT);
    AddGap(5.5);

    % 11 biot (10 BPM)
    t = 0:1/fs:21-1/fs;
    biot = 2.5*sin(2*pi*10/60*t);
    biot_rr = ones(1, length(t))*10;
    biot_a = ones(1, length(t))*2.5;
    biot = [biot zeros(1, round(19.25*fs))];
    biot_rr = [biot_rr zeros(1, round(19.25*fs))];
    biot_a = [biot_a zeros(1, round(19.25*fs))];
    biot = [biot, biot];
    biot_rr = [biot_rr, biot_rr];
    biot_a = [biot_a, biot_a];
    sig = [sig biot];
    rr = [rr biot_rr];
    a = [a biot_a];
    AddLabel(length(biot), Pattern.BIOT);
    AddGap(7.5);

    % 12 apnea (0 BPM)
    apnea = zeros(1, round(60*fs));
    sig = [sig, apnea]; 
    AddLabel(length(apnea), Pattern.APNEA);
    AddFeatures(length(apnea), 0, 0);
    
    if plot
        ShowSignal({sig, label}, fs);
    end
    
%     [p, q] = rat(fs_disp/fs);
%     sig = resample(sig, p, q);
%     label = ResampleLabels(label, fs_disp, fs);
    
    res.sig = sig';
    res.labels = label';
    res.source = 'reference';
    res.rr = rr';
    res.a_rr = a';
    res.var_rr = movvar(res.rr, 30*30);
    res.var_a_rr = movvar(res.a_rr, 30*30);

    function AddLabel(len, label_type)
        label_sec = zeros(1, len);
        label_sec(:) = label_type;

        label = [label label_sec];
    end

    function AddGap(len)
        gap = nan(1, round(len*fs));
        sig = [sig gap];
        rr = [rr gap];
        a = [a gap];
        AddLabel(length(gap), Pattern.UNDEFINED);   
    end

    function AddFeatures(len, rr_set, a_set)
        rr = [rr ones(1, len)*rr_set];
        a = [a ones(1, len)*a_set];
    end
end
