clc
close all

path(path,'utils/')
path(path,'utils/GraphCuts/')
path(path,'data/')

test_image = 'ICML';

switch test_image
    case 'Background'
        load ImgMask1
        c = 3;
        bg_index = 107; % frame index
        Img0 = ImgMask1{bg_index,1};
        mask=ImgMask1{bg_index, 2};
        bin=repmat(mask, [1, 1, 3]);
        Img0 = im2double(Img0);
        x=Img0.*bin;
        X0 = rgb2gray(x);
        x_gray = imresize(X0,[100 100],'bilinear');
        [imheight,imwidth] = size(x_gray);
    case 'Angio'
        Img = imread('../data/angiography.jpg');
        Img = im2double(Img); X = rgb2gray(Img);
        X = X(1:900,1:900); % make square
        X = 1-X; % invert pixels
        X = X-min(X(:)); X = X/max(X(:));
        Xsmall = X;
        Xsmall = imresize(Xsmall,[100 NaN],'bilinear');
        k = round(0.06*numel(Xsmall));
        c = 1;
        gamma = 0.04;
        supp = cluster_k_approx(Xsmall,k,c,gamma,false);
        x_gray = supp.*Xsmall;
        [imheight,imwidth] = size(x_gray);
    case 'ICML'
        Img = imread('../data/icml2.png');
        Img = im2double(Img);
        X = rgb2gray(Img);
        X = X - min(X(:)); X = X/max(X(:));
        Xsmall = imresize(X, [100 100]);
        c = 4;
        k = round(0.04*numel(Xsmall));
        gamma = 0.04;
        verbose = 0;
        supp = cluster_k_approx(Xsmall,k,c,gamma,verbose);
        x_gray = supp.*Xsmall;
        [imheight,imwidth] = size(x_gray); c = c-1;
end

% CoSaMP/ModelCS options
opts.selectAtom = 2;
opts.gamma = 0.1;
opts.tol=1e-4;
opts.maxIter=50;
opts.verbose=0;
opts.tol_early=0.005;
opts.debias = 0;
opts.hhs = 0;
opts.imheight = imheight;
opts.imwidth = imwidth;
opts.display_perf = 1;

n = numel(x_gray);
k = nnz(x_gray);
c = c+1;
%Mvec = 4;
Mvec = 2:0.2:7;
E = 50;
%E = 2;

error_spg = zeros(E,length(Mvec));
error_cosamp = zeros(E,length(Mvec));
error_graphcuts = zeros(E,length(Mvec));
error_clusters = zeros(E,length(Mvec));
error_structomp = zeros(E,length(Mvec));

time_spg = zeros(E,length(Mvec));
time_cosamp = zeros(E,length(Mvec));
time_graphcuts = zeros(E,length(Mvec));
time_clusters = zeros(E,length(Mvec));
time_structomp = zeros(E,length(Mvec));

trial_time_average = 0.0;

for ee=1:E
    fprintf('Trial %d\n', ee);
    trial_timer = tic;
    
    for mm = 1:length(Mvec)
        m = 2*round(Mvec(mm)*k/2);
        
        fprintf('  Oversampling factor %f (M = %d)\n', Mvec(mm), m);
        
        P = (1:n)';
        Q = randperm(n/2-1)+1;
        OMEGA = Q(1:m/2)';
        Af = @(z) A_f(z, OMEGA, P);
        At = @(z) At_f(z, n, OMEGA, P);
        
        sgvar = 0.025;
        y0 = Af(x_gray(:));
        yn = sgvar*randn(size(y0));
        y = y0 + yn;
                         
        % Signal reconstruction
        %%%% L1 optimization
        Afun = @(x,mode) Afun_spg(x,mode,Af,At);
        opts_spg = spgSetParms('optTol',1e-3,'iterations',opts.maxIter,'verbosity',0);
        disp('L1 minimization:')
        tstart = tic;
        xhat_spg = spg_bp(Afun, y, opts_spg);
        time_spg(ee,mm) = toc(tstart);
        error_spg(ee,mm) = norm(x_gray(:)-xhat_spg(:))^2/norm(x_gray(:))^2;
        snr_spg = -10*log10(error_spg(ee,mm));
        if opts.display_perf
            disp([snr_spg,time_spg(ee,mm)])
        end
        
        %%%% CoSaMP
        recflag = 'sparse';
        disp('CoSaMP:')
        tstart = tic;
        [xhat_cosamp,tcgs,~] = cosamp_cluster(recflag, y, k, c, Af, At, opts);
        time_cosamp(ee,mm) = toc(tstart);
        error_cosamp(ee,mm) = norm(x_gray(:)-xhat_cosamp(:))^2/norm(x_gray(:))^2;
        snr_cosamp = -10*log10(error_cosamp(ee,mm));
        if opts.display_perf
            disp([snr_cosamp,time_cosamp(ee,mm),tcgs])
        end
        
        %%%% GraphCuts
        %addpath ../../../StructOMP/Background/
        opts.SF = 0.25;
        opts.thr = 0.05;
        opts.iter = 20;
        opts.display_progress = 0;
        disp('Graph Cuts')
        tstart = tic;
        xhat_graphcuts = cosamp_graphcuts_fun(y,Af,At,k,opts);
        time_graphcuts(ee,mm) = toc(tstart);
        error_graphcuts(ee,mm) = norm(x_gray(:)-xhat_graphcuts(:))^2/norm(x_gray(:))^2;
        snr_graphcuts = -10*log10(error_graphcuts(ee,mm));
        if opts.display_perf
            disp([snr_graphcuts,time_graphcuts(ee,mm)])
        end
        
        %%% StructOMP
        g = c; K = k;
        H = opts.imheight; W = opts.imwidth;
        disp('Graph Sparsity:')
        Edge4=Edge4Index(H, W);
        [~, Bm]=GetBlocksMatrix(H,W,1);
        [~, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);
        mq=size(Bm, 1);
        lambda=1; cl0=2*lambda*g*log2(mq)+K;
        tstart = tic;
        [xhat_structomp,~,~,~] = GraphOMP_CS_chin(cl0,{Af,At},y, Bm, BCm, lambda, []);
        time_structomp(ee,mm) = toc(tstart);
        error_structomp(ee,mm) = norm(x_gray(:)-xhat_structomp(:))^2/norm(x_gray(:))^2;
        snr_structomp = -10*log10(error_structomp(ee,mm));
        if opts.display_perf
            disp([snr_structomp,time_structomp(ee,mm)])
        end
        
        %%% ModelCS
        recflag = 'cluster';
        tstart = tic;
        disp('Cluster model:')
        [xhat_clusters,tcgs,tmodel,tproj] = cosamp_cluster(recflag, y, k, c, Af, At, opts);
        time_clusters(ee,mm) = toc(tstart);
        error_clusters(ee,mm) = norm(x_gray(:)-xhat_clusters(:))^2/norm(x_gray(:))^2;
        snr_cluster = -10*log10(error_clusters(ee,mm));
        if opts.display_perf
            disp([snr_cluster,time_clusters(ee,mm),tcgs,tmodel,tproj])
        end
    end
    
    trial_time = toc(trial_timer);
    trial_time_average = (trial_time_average * (ee - 1) + trial_time) / ee;
    time_remaining = trial_time_average * (E - ee);
    fprintf('Trial %d took %f seconds. Projected time remaining: %f seconds\n', ee, trial_time, time_remaining);
end


tol = 0.05;

p_spg = mean((error_spg < tol),1);
p_cosamp = mean((error_cosamp < tol),1);
p_graphcuts = mean((error_graphcuts < tol),1);
p_structomp = mean((error_structomp < tol),1);
p_clusters = mean((error_clusters < tol),1);

figure(1), clf, hold on
box on
grid on
siz = 10;
axisfortex('','Oversampling ratio m/k','Probability of recovery')
plot(Mvec,p_spg,'k-x','LineWidth',2,'MarkerSize',siz)
plot(Mvec,p_cosamp,'r-o','LineWidth',2,'MarkerSize',siz)
plot(Mvec,p_graphcuts,'g-d','LineWidth',2,'MarkerSize',siz)
plot(Mvec,p_structomp,'m-d','LineWidth',2,'MarkerSize',siz)
plot(Mvec,p_clusters,'-*','LineWidth',2,'MarkerSize',siz)
%axis([2 6 0 1])
legend('SPGL1','CoSaMP','Graph-Cuts','StructOMP','ModelCS','Location','Best')

figure(2), clf, hold on
box on
grid on
siz = 10;
axisfortex('','Oversampling ratio m/k','Average recovery time')
plot(Mvec,(mean(time_spg,1)),'k-x','LineWidth',2,'MarkerSize',siz)
plot(Mvec,(mean(time_cosamp,1)),'r-o','LineWidth',2,'MarkerSize',siz)
plot(Mvec,(mean(time_graphcuts,1)),'g-d','LineWidth',2,'MarkerSize',siz)
plot(Mvec,(mean(time_structomp,1)),'m-d','LineWidth',2,'MarkerSize',siz)
plot(Mvec,(mean(time_clusters,1)),'-*','LineWidth',2,'MarkerSize',siz)
legend('SPGL1','CoSaMP','Graph-Cuts','StructOMP','ModelCS','Location','Best')


% file = fopen('mc_results.txt', 'w');
% fprintf(file, 'r p_exact_tree p_approx_tree p_cosamp p_l1\n');
% for mm = 1:length(Mvec)
%     fprintf(file, '%e %e %e %e %e\n', Mvec(mm), p_exact_tree(mm), p_approx_tree(mm), p_cosamp(mm), p_l1(mm));
% end
% fclose(file);
