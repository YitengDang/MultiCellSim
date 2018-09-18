function [period, t_onset] = periodicity_test_detailed(cells_hist,...
    period_ub, decimals)
    % for a confirmed periodic solution, checks within the last t_check+1
    % states when the first state revisit happened    
    
    % rounds the states up to a specified number of decimals (for finite
    % Hill system); if decimals = Inf, do not perform rounding
    if nargin<3
        decimals = Inf;
    end
    
    % uses period_ub as an upper bound for the period obtained from periodicity_test_short
    t_out = numel(cells_hist)-1;
    t_start = 2; % start checking at t=2
    for t1=t_start:t_out
        %disp(t1);
        if decimals<Inf
            cells_1 = round(cells_hist{t1+1}, decimals);
        else
            cells_1 = cells_hist{t1+1};
        end
        t_start2 = max(t1-period_ub, 0); % do not start before t=0
        t_out2 = t1-2; % min. period = 2, so terminate 2 steps before end
        for t=t_start2:t_out2 
            %disp(t);
            if decimals<Inf
                cells_2 = round(cells_hist{t+1}, decimals);
            else
                cells_2 = cells_hist{t+1};
            end
            if all(all(cells_2==cells_1))
                period = t1-t;
                t_onset = t;
                %fprintf('t1=%d, period %d \n', t_onset, period);
                return
            end
        end
    end
    period = Inf;
    t_onset = Inf;
    warning('No period found in periodicity_test_detailed!')
    disp('Error: no period found');
end