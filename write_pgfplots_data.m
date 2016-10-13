function write_pgfplots_data(filename, header, data)
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', header);
    fclose(fid);
    dlmwrite(filename, data, 'delimiter', ' ', 'precision', '%.6e', '-append');
end

