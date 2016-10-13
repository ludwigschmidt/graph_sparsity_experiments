function [xhat,tcgs,tmodel,tproj] = cosamp_cluster(recflag, b, k, c, A, At, opts)

%%%%%%%%%%%

temp = At(b);
n = length(temp);
err = 1; 
aa= zeros(n,1); 
s_cosamp = zeros(n,1);
alph = opts.selectAtom;

iterCnt = 0; tcgs = 0; tmodel = 0; tproj = 0;

while ((err>opts.tol)&&(iterCnt<=opts.maxIter))
    iterCnt = iterCnt + 1;
    
    %disp(iterCnt);
    
    resid = b - A(s_cosamp);
    proxy = At(resid); 
    %------Estimation
    switch recflag
        case 'sparse'
            supp = thresh(proxy,round(alph*k));
        case 'cluster'
            gamma = opts.gamma;
            Xgrid = reshape(proxy,[],opts.imwidth); Xgrid = Xgrid.^2;
            verbose = opts.verbose;
            tstart_model = tic;
            [supp,t] = cluster_k_approx(Xgrid,round(alph*k),round(alph*c),gamma,verbose);
            tmodel = tmodel + toc(tstart_model); tproj = tproj + t;
            supp = supp(:);
    end
    
    if opts.hhs
        tt = find(supp);
        samples = resid;
    else
        tt = union(find(ne(s_cosamp,0)),find(supp));
        samples = b;
    end
    
    %------Least-squares
    PP_tt = @(z) A_I(A,z,tt,n);
    PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
    qq = PP_transpose_tt(samples);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    tstart_cgs = tic;
    w = cgsolve(PPtranspose_PP_tt,qq, 1e-4, 20, 0);
    tcgs = tcgs + toc(tstart_cgs);
    %w = cgs(PPtranspose_PP_tt,qq, 1e-4, 100);
    bb= 0*s_cosamp; bb(tt)= w;
    
    %------OPTIONAL: Merge
    if opts.hhs
        bb = s_cosamp + bb;
    end
    
    %------Prune
    switch recflag
        case 'sparse'
            supp = thresh(bb,k);
        case 'cluster'
            gamma = opts.gamma;
            alph = opts.selectAtom;
            Xgrid = reshape(bb,[],opts.imwidth);Xgrid = Xgrid.^2;
            verbose = opts.verbose;
            tstart_model = tic;
            [supp,t] = cluster_k_approx(Xgrid,k,c,gamma,false);
            tmodel =  tmodel + toc(tstart_model); tproj = tproj + t;
            supp = supp(:);
    end
    bb=bb(:);
    s_cosamp = supp.*bb;
    
    aa = s_cosamp;
    
    errTemp = norm(A(s_cosamp)-b)/norm(b);
    if (opts.verbose)
        fprintf('Iter: %d Err: %f\n',iterCnt,errTemp);
    end
    
    if (abs((err-errTemp)/err) <= opts.tol_early || errTemp <= opts.tol_early) %Early termination condition
       % err = errTemp;
        if (opts.verbose)
            fprintf('Terminating.\n');
        end
        break 
    end
    err = errTemp;
end

if opts.debias
    tt = find(supp);
    PP_tt = @(z) A_I(A,z,tt,n);
    PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
    qq = PP_transpose_tt(b);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    w = cgsolve(PPtranspose_PP_tt,qq, 1e-4, 100, 0);
    xhat = 0*s_cosamp; xhat(tt)= w;
else
    xhat = aa;
end

end

function supp = thresh(x,k)
[~,ww] = sort(abs(x(:)),'descend');
supp = 0*x;
supp(ww(1:k)) = 1;
end