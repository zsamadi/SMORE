function writeH5X(filename,data, indices, indptr)
% % Define vectors
% data = ones(100, 1);    % Replace with your actual data vector
% indices = ones(100, 1); % Replace with your actual indices vector
% indptr = ones(100, 1);   % Replace with your actual indptr vector

if isfile(filename)
    delete(filename)
end

% Create group 'X' and write datasets into it
h5create(filename, '/X/data', size(data));
h5write(filename, '/X/data', data);

h5create(filename, '/X/indices', size(indices));
h5write(filename, '/X/indices', indices);

h5create(filename, '/X/indptr', size(indptr));
h5write(filename, '/X/indptr', indptr);

disp('Data written successfully to HDF5 file.');