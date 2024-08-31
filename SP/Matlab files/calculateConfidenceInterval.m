function confidenceInterval = calculateConfidenceInterval(data)
    % Number of bootstrap resamples
    numResamples = 10000;
    
    % Confidence level
    alpha = 0.95;

    % Initialize an array to store resampled means
    resampledMeans = zeros(numResamples, 1);

    for i = 1:numResamples
        % Generate a resample by sampling with replacement
        resample = datasample(data, length(data));

        % Compute the mean of the resample
        resampledMeans(i) = mean(resample);
    end

    % Sort the resampled means
    sortedMeans = sort(resampledMeans);

    % Calculate the lower and upper percentiles for the confidence interval
    lowerPercentile = (1 - alpha) / 2 * 100;
    upperPercentile = 100 - lowerPercentile;

    % Calculate the confidence interval
    confidenceInterval = [sortedMeans(round(lowerPercentile/100 * numResamples)), sortedMeans(round(upperPercentile/100 * numResamples))];
end
