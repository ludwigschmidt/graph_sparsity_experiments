function [supp,tproj] = cluster_k_approx(coeffs,k,c,gamma,verbose)


% Inputs: 
% coeffs: signal coefficients in grid-format
% k: target sparsity
% c: target #clusters
% verbose: display progress

% Output
% supp: optimal cluster-sparse support
% Note that supp can have between [k (1+gamma)*k] nonzeros

use_cpp = true;

if use_cpp
    opts = [];
    opts.pruning = 'gw';
    opts.max_num_iter = 10;
    if verbose
        opts.verbose = 1;
    end
    tproj_start = tic;
    [supp, ~] = cluster_grid_pcst_binsearch(coeffs, c, k, round((1 + gamma) * k), opts);
    tproj = toc(tproj_start);
else
    % perform binary search on lambda to get a good k
    lambda_lo = 0; 
    lambda_hi = max(coeffs(:)); 
    %coeffs_sort = sort(coeffs(:),'descend');
    %lambda_hi = coeffs_sort(2*k);
    lambda_tol = eps;

    ii = 0; tproj = 0;

    while lambda_hi-lambda_lo > lambda_tol
        ii=ii+1;

        lambda = 0.5*(lambda_hi + lambda_lo);

        tproj_start = tic;
        [supp,~] = cluster_grid_pcst(coeffs,c,lambda);
        ksupp = nnz(supp);

        tproj = tproj+toc(tproj_start);

        if ksupp >= k && ksupp <= round((1+gamma)*k)
            %disp([ksupp >= k, ksupp <= round((1+gamma)*k)])
            %disp('Ending bin search')
            break;
        end
        if ksupp > ((1+gamma)*k) 
            lambda_lo = lambda;
        else
            lambda_hi = lambda;
        end

        if verbose
            fprintf('Sparsity = %d, Lambda_lo = %f, Lambda_hi = %f\n',ksupp,lambda_lo,lambda_hi);
        end

    end
end

end