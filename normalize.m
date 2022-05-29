function err_norm  = normalize(errors)

max_viol = 0.001 + max(errors);      % Max violation
[err_size,~] = size(errors);         % Error matrix size
max_viol_matrix = repmat(max_viol,err_size,1);  % Max violation matrix
err_norm1 = errors./max_viol_matrix;        % Normalized errors
err_norm = sum(err_norm1,2);               % Total sum
end