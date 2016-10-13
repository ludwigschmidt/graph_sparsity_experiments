function write_phasetransition_data(filename, oversamplings, p_models, p_structomp, p_graphcuts, p_cosamp, p_spgl1)
    num_rows = numel(oversamplings);
    data = zeros(num_rows, 6);
    data(:,1) = oversamplings;
    data(:,2) = p_models;
    data(:,3) = p_structomp;
    data(:,4) = p_graphcuts;
    data(:,5) = p_cosamp;
    data(:,6) = p_spgl1;
    write_pgfplots_data(filename, 'oversampling p_models p_structomp p_graphcuts p_cosamp p_spgl1', data);
end

