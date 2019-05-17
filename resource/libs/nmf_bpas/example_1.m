% ------------------------------ Usage examples ------------------------------
m = 500;
n = 100;
k = 100;

A = rand(m,n);

% Uncomment only one -------------------------------------------------------------------------
tic
[W_nmf,H_nmf,iter_nmf,HIS_nmf]=nmf(A,k);
tnmf = toc();
% [W,H,iter,HIS]=nmf(A,k,'verbose',2);

tic
[W_as,H_as,iter_as,HIS_as]=nmf(A,k,'verbose',1,'nnls_solver','as');
t_as = toc();

tic
[W_bp,H_bp,iter_bp,HIS_bp]=nmf(A,k,'verbose',1,'nnls_solver','bp');
t_bp = toc();

% [W,H,iter,HIS]=nmf(A,k,'verbose',1,'type','sparse');
% [W,H,iter,HIS]=nmf(A,k,'verbose',1,'type','sparse','nnls_solver','bp','alpha',1.1,'beta',1.3);
% [W,H,iter,HIS]=nmf(A,k,'verbose',2,'type','plain','w_init',rand(m,k));

