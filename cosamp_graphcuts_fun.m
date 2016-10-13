function xhat= cosamp_graphcuts_fun(y, Af, At, s, opts)
%---
y = y(:); % 
imwidth = opts.imwidth;
imheight = opts.imheight;
%thr = opts.thr;
N = imheight*imwidth;
%---
%---
aa = zeros(N,opts.iter); % stores current sparse estimate
kk=1; % current MP iteration
maxiter_cg= 500;
verbose= 0;
tol= 1e-3;
tic
while le(kk,opts.iter),
    rr = y - Af(aa(:,kk));
    proxy = At(rr);
    %---
    if kk>1
        yyi=proxy;
        yyi = yyi+aa(:,kk);
        yyisort = sort(abs(yyi),'descend'); 
        thr = yyisort(2*s);
        
        L1= zeros(size(yyi));
        L2= zeros(size(yyi));
        L1(abs(yyi)>=thr)= 1; 
        L1(abs(yyi)<thr)= 0;
        L2(abs(yyi)<thr)= 1;
        L2(abs(yyi)>=thr)= 0;
        Dc = zeros(imheight,imwidth,2); 
        % data fidelity cost (size(Dc)=imheight*imwidth*L) 
        Dc(:,:,1) = reshape(L1,imheight,imwidth);
        Dc(:,:,2) = reshape(L2,imheight,imwidth);
        % smoothness cost (size(SC)=L*L)
        Sc = ones(2) - eye(2);   
        % Graph cut 
        gch = GraphCut('open', Dc, opts.SF*Sc);
        [gch L] = GraphCut('expand',gch);
        gch = GraphCut('close', gch);       
        if opts.display_progress
            figure(10), clf
            imshow(L,[])
            pause(0.3), shg
        end
        % Find sparsity pattern after enforcing graphical model.
        tt= find(L(:)==1);
    else
        % used only for the first MP iteration
        [~,ww]= sort(abs(proxy),'descend');
        tt= union(find(ne(aa(:,kk),0)),ww(1:(2*s)));
    end

    PP_tt = @(z) A_I(Af,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
    qq = PP_transpose_tt(y);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    x = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter_cg, verbose);

    bb= zeros(N,1);
    bb(tt)= x;
    %---
    [~,ww2]= sort(abs(bb),'descend');
    %---
    kk= kk+1;
    aa(ww2(1:s),kk)= bb(ww2(1:s)); 
end
xhat= aa(:,end);