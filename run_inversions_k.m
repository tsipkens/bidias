
% RUN_INVERSIONS_K  Time Bayesian and Twomey-Markowski optimization routines. 
% Author: Timothy Sipkens, 2020-02-22
%=========================================================================%

clear t;

[x_two_mh,Sf_two_mh,out_two_mh] = ...
    optimize.twomark_op(A,b,Lb,grid_x,...
    x_init,35,[1e2,1e-5],x0,'Buckley');
eps.two_mh = norm(x0-x_two_mh);

%%
nt = 25;
disp('Performing time repeats...');
tools.textbar(0);
for tt=1:nt
    tic;
    x_two_mh = invert.twomark(A,b,Lb,grid_x,...
        x_init,35,'Buckley',1/Sf_two_mh);
    t(tt).two_mh = toc;
    disp('Completed Twomey-Markowski.');
    
    eps.two_mh = norm(x0-x_two_mh);
    
    
    tic;
    invert.em(Lb*A,Lb*b,x_init,3,x0);
    t(tt).em = toc;
    
    
    tic;
    invert.tikhonov(...
        Lb*A,Lb*b,lambda_tk1,1,grid_x,[]);
    t(tt).tk1 = toc;
    
    
    
    tic;
    invert.exp_dist(...
        Lb*A,Lb*b,lambda_ed_lam,Gd,...
        grid_x);
    t(tt).ed = toc;
    
    disp('Performing time repeats...');
    tools.textbar(0);
    tools.textbar(tt/nt);
end

tmean.two_mh = mean([t.two_mh]);
tmean.em = mean([t.em])./tmean.two_mh;
tmean.tk1 = mean([t.tk1])./tmean.two_mh;
tmean.ed = mean([t.ed])./tmean.two_mh;
tmean.two_mh = 1;

