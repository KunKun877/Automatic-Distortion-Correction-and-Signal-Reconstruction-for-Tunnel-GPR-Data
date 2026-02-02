function consistency = fun_trace_detect(field)
    data = field;
    win_size = 1;
    [n_samples, n_traces] = size(data); % Correctly get dimensions
    
    consistency = ones(n_traces, 1);
    
    loop_start_idx = win_size + 1;
    loop_end_idx = n_traces - win_size;
    
    % Calculate neighborhood consistency for each trace
    for i = loop_start_idx:loop_end_idx
%         start_idx = i;
%         end_idx = i+1;
        start_idx = max(1, i - win_size);
        end_idx = min(n_traces, i + win_size);
    
    
        neighborhood_traces = data(:, start_idx:end_idx);
        
%         ref_trace = median(neighborhood_traces, 2); % Median along dimension 2 (across columns/traces)
        ref_trace = mean(neighborhood_traces, 2);
        
        current_trace = data(:, i);
        consistency(i) = corr(current_trace, ref_trace);

end