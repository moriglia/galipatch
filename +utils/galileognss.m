function [f_c, B] = galileognss(signal_name)
    %% Galileo Frequencies and Bandwidths
    carrier_frequency = containers.Map;
    carrier_frequency('E1')     = 1575.420e6;
    carrier_frequency('E5a')    = 1176.450e6;
    carrier_frequency('E5b')    = 1207.150e6;
    carrier_frequency('E6')     = 1278.750e6;

    signal_bandwidth = containers.Map;
    signal_bandwidth('E1')      = 24.552e6;
    signal_bandwidth('E5a')     = 20.460e6;
    signal_bandwidth('E5b')     = 20.460e6;
    signal_bandwidth('E6')      = 40.920e6;

    % f0 = [1575.420 1176.450 1207.150 1278.750]*1e6;
    % B = [ 24.552 20.460 20.460 40.920]*1e6;
    
    if signal_name == "all"
        f_c = carrier_frequency;
        B = signal_bandwidth;
    else
        f_c = carrier_frequency(signal_name);
        B = signal_bandwidth(signal_name);
    end
end