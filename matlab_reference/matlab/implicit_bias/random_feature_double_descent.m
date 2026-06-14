clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('random_feature_dd.fig');
catch
    disp('missing figure file')
end


d = 5;
ntrain = 200;
ntest = 200;
m = 100;
sigma = 2;
w0 = randn(d,1);
Xtrain = randn(ntrain,d);
Xtrain = Xtrain ./ sqrt(repmat(sum(Xtrain.^2,2),1,d));
ytrain = 1./( .25 + (Xtrain*w0).^2) + sigma * randn(ntrain,1);

Xtest = randn(ntest,d);
Xtest = Xtest ./ sqrt(repmat(sum(Xtest.^2,2),1,d));
ytest = 1./( .25 + (Xtest*w0).^2) + sigma * randn(ntest,1);

W_rffs = randn(m,d);
W_rffs = W_rffs ./ sqrt(repmat(sum(W_rffs.^2,2),1,d));
Phi_rffs_train = max(Xtrain * W_rffs',0);
Phi_rffs_test = max(Xtest * W_rffs',0);


V = ( Phi_rffs_train' * Phi_rffs_train ) \ ( Phi_rffs_train' * ytrain);
ytest_pred = Phi_rffs_test * V;
plot(ytest,ytest_pred,'x');

ms = 0:1:800;

W_rffs = randn(max(ms),d);
W_rffs = W_rffs ./ sqrt(repmat(sum(W_rffs.^2,2),1,d));


% regularization by feature cardinality

for im = 1:length(ms);
    m = ms(im);
    Phi_rffs_train = max(Xtrain * W_rffs(1:m,:)',0);
    Phi_rffs_test = max(Xtest * W_rffs(1:m,:)',0);
    
    
    V = ( Phi_rffs_train' * Phi_rffs_train + 1e-12 * eye(m) * ntrain) \ ( Phi_rffs_train' * ytrain);
    ytest_pred = Phi_rffs_test * V;
    ytrain_pred = Phi_rffs_train * V;
    train_errors(im) = mean( (ytrain - ytrain_pred).^2);
    test_errors(im) = mean( (ytest - ytest_pred).^2);
end


% regularization by square norm

lambdas = 10.^[-12:.25:4];
m = max(ms);
Phi_rffs_train = max(Xtrain * W_rffs(1:m,:)',0);
Phi_rffs_test = max(Xtest * W_rffs(1:m,:)',0);

for ilambda = 1:length(lambdas);
    lambda = lambdas(ilambda);
    V = ( Phi_rffs_train' * Phi_rffs_train + lambda * eye(m) * ntrain) \ ( Phi_rffs_train' * ytrain);
    ytest_pred = Phi_rffs_test * V;
    ytrain_pred = Phi_rffs_train * V;
    train_errors_lambda(ilambda) = mean( (ytrain - ytrain_pred).^2);
    test_errors_lambda(ilambda) = mean( (ytest - ytest_pred).^2);
end



subplot(1,2,1)
plot(1:ntrain,train_errors(1:ntrain),'b','linewidth',2); hold on; plot(1:ntrain,test_errors(1:ntrain),'r','linewidth',2); hold off
set(gca,'fontsize',20);
legend('train','test','location','northwest');
axis([0 ntrain-1 0 19])
xlabel('number of random features')
ylabel('test error')
title('Feature cardinality','FontWeight','normal');

subplot(1,2,2)
plot(-log10(lambdas),train_errors_lambda,'b','linewidth',2); hold on; plot(-log10(lambdas),test_errors_lambda,'r','linewidth',2); hold off
set(gca,'fontsize',20);
legend('train','test','location','northwest');
axis([min(-log10(lambdas)) max(-log10(lambdas)) 0 19])
xlabel('-log_{10}(\lambda)')
ylabel('test error')
title('Regularization parameter','FontWeight','normal');


try
    print('-depsc', 'ushaped_rf.eps');
        close(ccc)
catch
    disp('missing figure file')
end


try
    ccc=openfig('doubledescent_top.fig');
catch
    disp('missing figure file')
end


plot(ms,train_errors,'b','linewidth',2); hold on; plot(ms,test_errors,'r','linewidth',2); hold off
set(gca,'fontsize',20);
legend('train','test','location','southeast');
axis([0 max(ms) 0 19])
xlabel('number of random features')
ylabel('test error')
title('Feature cardinality','FontWeight','normal');




try
    print('-depsc', 'doubledescent_top.eps');
    close(ccc)
catch
    disp('missing figure file')
end

