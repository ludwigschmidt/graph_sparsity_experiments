function [xest, input, tprune, tcgs] = GraphOMP_CS_chin(cl0,A,y, Bm, BC, lamada, input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Graph structured sparsity
%%% Input:
%%%   cl0:  the coding complexity
%%%   Phi:  projection matrix
%%%   y:   measurements  
%%%   Bm:  the defined block matrix, where the i-th row denotes all index
%%%        entires in sparse x included in the i-th block.
%%%   BC:  the connected relations between blocks, where the i-th row
%%%        denotes all index of blocks connected to the i-th blocks 
%%%   lambda: weights in the coding complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output:
%%%   xest: the unknow structured sparse data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with
%%%%   Structured Sparsity", Rutgers University.
%%%%   By Junzhou Huang, jzhuang@cs.rutgers.edu
%%%%   Jan 2009, Updated Dec 20, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modified by Chin to work with function handles

Af = A{1};
At = A{2};
m = length(y);
n = length(At(y));

tprune = 0; tcgs = 0;

if isempty(input)
    y_r = y;
    input.supp = [];input.Blist=[];
    input.cl=0; input.RegionNum=0; 
    input.RegionIndex=[];
    x_p=[];
else
    y_r=y - Af(input.x);
    x_p=input.x_p;
end

maxiter_lsqr=100;

in = 0;
while in<10000
   in = in+1;
   cv = abs( At(y_r) ); cv = cv';
   tstart_prune = tic;
   [input.supp, input.Blist, input.cl, input.RegionNum,input.RegionIndex, newid] = GraphPruneMex(cv, Bm, BC, lamada,input.supp, input.Blist, input.cl, input.RegionNum, input.RegionIndex);
   tprune = tprune + toc(tstart_prune);
   curr_index=input.supp;
%    
%    Phi_x = Phi(:,curr_index);
%    x0=zeros(length(curr_index),1);
%    idd=~ismembc(curr_index, newid);
%    indx0=find(idd==1);
%    x0(indx0')=x_p;   
%    [x_p, flag] = lsqr_gp(Phi_x, y, curr_index', OptTol, maxiter_lsqr, [], [], x0);   
%    y_r = y - Phi_x*x_p;
%    
   PP_tt = @(z) A_I(Af,z,curr_index,n);
   PP_transpose_tt = @(z) A_I_transpose(At,z,curr_index);
   qq = PP_transpose_tt(y);
   PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
   tstart_cgs = tic;
   x_p = cgsolve(PPtranspose_PP_tt,qq, 1e-6, maxiter_lsqr, 0);
   tcgs = tcgs + toc(tstart_cgs);
   y_r = y - PP_tt(x_p);
   
   %norm_save(in) = norm(y_r);
%    if mod(in,100) == 0
%        disp(in);
%    end
   if input.cl>cl0
       break;
   end
end
xest=zeros(n, 1);
xest(curr_index,1)=x_p;
input.x=xest;
input.x_p=x_p;
return