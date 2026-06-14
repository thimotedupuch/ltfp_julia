clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('cg_lasso.fig');
catch
    disp('missing figure file')
end

d = 1000;
n = 100;
ntest = 1000;
k = 5;
sigma_noise = .5;
 
X = randn(n,d);
Xtest = randn(ntest,d);
beta_ast = zeros(d,1);
beta_ast(1:k) = sign(randn(k,1));



% sparse

y = X * beta_ast + sigma_noise * randn(n,1);
ytest = Xtest * beta_ast + sigma_noise * randn(ntest,1);



maxiter  = 50;
beta = zeros(d,1);
L = max(eig(X'*X/n));
steps_size_prox = 1/L;
Ls = diag(X'*X/n) ;
u = zeros(n,1);

for iter=1:maxiter
    
    % conditional gradient
    betas(:,iter) = beta;
    
    perftrain(iter) = mean( (y-X*beta).^2);
    perftest(iter) = mean( (ytest-Xtest*beta).^2);
    
    % boosting
    grad = u - y;
    [a,b] = max( abs( X' * grad ) );
    betaloc = zeros(d,1);
    if  ( X(:,b)'* grad  > 0 ) , betaloc(b) = -1; else betaloc(b) = 1; end
    uloc = X * betaloc;
    alpha = - ( uloc'* grad ) / ( uloc' * uloc );
    u = u + uloc * alpha;
    beta = beta + betaloc * alpha;    
    
end



subplot(2,2,1);
plot(0:maxiter-1,log(perftrain),'b','linewidth',2); hold on;
plot(0:maxiter-1,log(perftest),'r','linewidth',2); hold off
set(gca,'fontsize',24)
legend('train','test');
xlabel('iterations');
ylabel('log(errors)');
title('Errors - sparse','FontWeight','normal')
axis([0 maxiter-1  -5  3]);

  

subplot(2,2,2);
plot(0:maxiter-1,betas','linewidth',2);
xlabel('iterations')
set(gca,'fontsize',24)
ylabel('weights');
axis([0 maxiter-1  -1.2 1.2]);
title('Weights - sparse','FontWeight','normal')

 
temp = randn(d,d);
[u,e] = eig( temp+temp');
beta_ast = u * beta_ast;

y = X * beta_ast + sigma_noise * randn(n,1);
ytest = Xtest * beta_ast + sigma_noise * randn(ntest,1);



% non - sparse
maxiter  = 50;
beta = zeros(d,1);
L = max(eig(X'*X/n));
steps_size_prox = 1/L;
Ls = diag(X'*X/n) ;
u = zeros(n,1);

for iter=1:maxiter
    
    % conditional gradient
    betas(:,iter) = beta;
    
    perftrain(iter) = mean( (y-X*beta).^2);
    perftest(iter) = mean( (ytest-Xtest*beta).^2);
    
    % boosting
    grad = u - y;
    [a,b] = max( abs( X' * grad ) );
    betaloc = zeros(d,1);
    if  ( X(:,b)'* grad  > 0 ) , betaloc(b) = -1; else betaloc(b) = 1; end
    uloc = X * betaloc;
    alpha = - ( uloc'* grad ) / ( uloc' * uloc );
    u = u + uloc * alpha;
    beta = beta + betaloc * alpha;    
    
end

subplot(2,2,3);
plot(0:maxiter-1,log(perftrain),'b','linewidth',2); hold on;
plot(0:maxiter-1,log(perftest),'r','linewidth',2); hold off
set(gca,'fontsize',24)
legend('train','test');
xlabel('iterations');
ylabel('log(errors)');
title('Errors - non sparse','FontWeight','normal')
axis([0 maxiter-1  -5  3]);



 

subplot(2,2,4);
plot(0:maxiter-1,betas','linewidth',2);
xlabel('iterations')
set(gca,'fontsize',24)
ylabel('weights');
axis([0 maxiter-1  -1.2 1.2]);
title('Weights - non sparse','FontWeight','normal')



try
    print('-depsc', 'cg_lasso.eps');
    close(ccc)
catch
    disp('missing figure file')
end


