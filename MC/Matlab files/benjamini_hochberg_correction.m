function p_values_adj = benjamini_hochberg_correction(p_values)
    % Calculate ranks and p_m_over_k
    ranks = tiedrank(p_values, 'last');
    p_m_over_k = p_values * length(p_values) ./ ranks;

    % Initialize the adjusted p-values
    p_values_adj = zeros(1, length(p_values));

    % Calculate adjusted p-values
    for i = 1:length(p_values)
        % Find the rank
        tmp_rank = ranks(i);

        % Get all the p_m_over_k that are greater or equal to this rank
        % and get the minimum value
        p_values_adj(i) = min(1, min(p_m_over_k(ranks >= tmp_rank)));
    end
end
