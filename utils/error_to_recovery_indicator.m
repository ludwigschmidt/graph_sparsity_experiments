function errors = error_to_recovery_indicator(errors, threshold)
    recovered = errors <= threshold;
    errors(:) = 0;
    errors(recovered) = 1;
end

