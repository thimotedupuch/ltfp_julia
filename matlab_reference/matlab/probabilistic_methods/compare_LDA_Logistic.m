clear all
seed = 2;
rand('state',seed);
randn('state',seed);

try
    ccc=openfig('compare_LDA_Logistic.fig');
catch
    disp('missing figure file')
end


ts = 0:.1:1;
nrep = 100;
d = 50;

for irep = 1:nrep
    % load data
    n1 = 50;
    n0 = 50;
    std_noise = 1;
    mu1ast = randn(d,1)/sqrt(d)*2;
    mu0ast = randn(d,1)/sqrt(d)*2;
    noise1_train = randn(n1,d);
    noise0_train = randn(n0,d);
    noise1_test = randn(n1,d);
    noise0_test = randn(n0,d);
    
    for it=1:length(ts)
        Sigma_common_sqrt1 = expm(-ts(it)*randn(d,d)/8);
        Sigma_common_sqrt0 = expm(-ts(it)*randn(d,d)/8);
        
        X1 = ones(n1,1) * mu1ast' + std_noise * noise1_train * Sigma_common_sqrt1;
        X0 = ones(n0,1) * mu0ast' + std_noise * noise0_train * Sigma_common_sqrt0;
        X = [ X1; X0];
        y = [ ones(n1,1); zeros(n0,1) ];
        n = n1 + n0;
        
        X1_test = ones(n1,1) * mu1ast' + std_noise *  noise1_test * Sigma_common_sqrt1;
        X0_test = ones(n0,1) * mu0ast' + std_noise * noise0_test * Sigma_common_sqrt0;
        X_test = [ X1_test; X0_test];
        y_test = [ ones(n1,1); zeros(n0,1) ];
        
        
        % LDA
        % maximum likelihood estimation
        mu1 = sum( X1, 1)'/n1;
        mu0 = sum( X0, 1)'/n0;
        
        Sigma1 = ( X1 - ones(n1,1)*mu1' )' *  ( X1 - ones(n1,1)*mu1' ) / n1;
        Sigma0 = ( X0 - ones(n0,1)*mu0' )' *  ( X0 - ones(n0,1)*mu0' ) / n0;
        
        Sigma = n0/n * Sigma0 + n1/n * Sigma1;
        
        % computing discriminative affine function
        w_LDA  = Sigma \ ( mu1 - mu0 );
        b_LDA = - .5 * ( mu1 + mu0 )' * w_LDA + log ( n1 / n0 );
        
        
        % compute error rates
        error_train_LDA(irep,it) = mean( abs( ( sign( X * w_LDA + b_LDA ) + 1 ) / 2 - y ) )
        error_test_LDA(irep,it) = mean( abs( ( sign( X_test * w_LDA + b_LDA ) + 1 ) / 2 - y_test ) )
        
        
        % LDA - diagonal
        % maximum likelihood estimation
        mu1 = sum( X1, 1)'/n1;
        mu0 = sum( X0, 1)'/n0;
        
        Sigma1 = ( X1 - ones(n1,1)*mu1' )' *  ( X1 - ones(n1,1)*mu1' ) / n1;
        Sigma0 = ( X0 - ones(n0,1)*mu0' )' *  ( X0 - ones(n0,1)*mu0' ) / n0;
        
        Sigma = diag(diag(n0/n * Sigma0 + n1/n * Sigma1));
        
        % computing discriminative affine function
        w_LDA  = Sigma \ ( mu1 - mu0 );
        b_LDA = - .5 * ( mu1 + mu0 )' * w_LDA + log ( n1 / n0 );
        
        
        % compute error rates
        error_train_LDA_diag(irep,it) = mean( abs( ( sign( X * w_LDA + b_LDA ) + 1 ) / 2 - y ) )
        error_test_LDA_diag(irep,it) = mean( abs( ( sign( X_test * w_LDA + b_LDA ) + 1 ) / 2 - y_test ) )
        
        
        
        % LOGISTIC REGRESSION
        Xint = [ X ones(n,1) ];   % add constant variable to x
        theta = zeros(d+1,1);
        normgradient = Inf;
        k=1;
        while (normgradient > 1e-10 & k < 100)
            eta = Xint * theta;
            mu = 1./(1+exp(-eta));
            gradient = Xint'*(y-mu) - 1e-12 * theta;
            H = - Xint' * ( repmat( mu.*(1-mu) , 1 , d+1 ) .* Xint ) - 1e-12 * eye(d+1);
            theta = theta - H \ gradient;
            if any(isnan(theta)), pause; end
            normgradient = norm(gradient);
            fprintf('k = %d - norm of gradient = %e\n',k,normgradient);
            k=k+1;
        end
        theta_LOGISTIC = theta;
        w_LOGISTIC = theta(1:d);
        b_LOGISTIC = theta(d+1);
        % compute error rates
        error_train_LOGISTIC(irep,it) = mean( abs( ( sign( X * w_LOGISTIC + b_LOGISTIC ) + 1 ) / 2 - y ) )
        error_test_LOGISTIC(irep,it) = mean( abs( ( sign( X_test * w_LOGISTIC + b_LOGISTIC ) + 1 ) / 2 - y_test ) )
        
    end
end
subplot(1,2,1);
plot(ts,mean(error_test_LDA_diag,1),'r','linewidth',2); hold on;
plot(ts,mean(error_test_LDA,1),'b','linewidth',2); hold on;
plot(ts,mean(error_test_LOGISTIC,1),'k','linewidth',2); hold off;
legend('LDA Diag','LDA','logistic','location','southeast');
set(gca,'fontsize',16)
xlabel('lack of independence');
ylabel('error rate')
axis([0 1 0 .23]);
title('d = 50, n = 50', 'fontweight','normal');




seed = 2;
rand('state',seed);
randn('state',seed);

ts = 0:.1:1;
nrep = 50;
d = 50;

for irep = 1:nrep
    % load data
    n1 = 500;
    n0 = 500;
    std_noise = 1;
    mu1ast = randn(d,1)/sqrt(d)*2;
    mu0ast = randn(d,1)/sqrt(d)*2;
    noise1_train = randn(n1,d);
    noise0_train = randn(n0,d);
    noise1_test = randn(n1,d);
    noise0_test = randn(n0,d);
    
    for it=1:length(ts)
        Sigma_common_sqrt1 = expm(-ts(it)*randn(d,d)/8);
        Sigma_common_sqrt0 = expm(-ts(it)*randn(d,d)/8);
        
        X1 = ones(n1,1) * mu1ast' + std_noise * noise1_train * Sigma_common_sqrt1;
        X0 = ones(n0,1) * mu0ast' + std_noise * noise0_train * Sigma_common_sqrt0;
        X = [ X1; X0];
        y = [ ones(n1,1); zeros(n0,1) ];
        n = n1 + n0;
        
        X1_test = ones(n1,1) * mu1ast' + std_noise *  noise1_test * Sigma_common_sqrt1;
        X0_test = ones(n0,1) * mu0ast' + std_noise * noise0_test * Sigma_common_sqrt0;
        X_test = [ X1_test; X0_test];
        y_test = [ ones(n1,1); zeros(n0,1) ];
        
        
        % LDA
        % maximum likelihood estimation
        mu1 = sum( X1, 1)'/n1;
        mu0 = sum( X0, 1)'/n0;
        
        Sigma1 = ( X1 - ones(n1,1)*mu1' )' *  ( X1 - ones(n1,1)*mu1' ) / n1;
        Sigma0 = ( X0 - ones(n0,1)*mu0' )' *  ( X0 - ones(n0,1)*mu0' ) / n0;
        
        Sigma = n0/n * Sigma0 + n1/n * Sigma1;
        
        % computing discriminative affine function
        w_LDA  = Sigma \ ( mu1 - mu0 );
        b_LDA = - .5 * ( mu1 + mu0 )' * w_LDA + log ( n1 / n0 );
        
        
        % compute error rates
        error_train_LDA(irep,it) = mean( abs( ( sign( X * w_LDA + b_LDA ) + 1 ) / 2 - y ) )
        error_test_LDA(irep,it) = mean( abs( ( sign( X_test * w_LDA + b_LDA ) + 1 ) / 2 - y_test ) )
        
        
        
        % LDA - diagonal
        % maximum likelihood estimation
        mu1 = sum( X1, 1)'/n1;
        mu0 = sum( X0, 1)'/n0;
        
        Sigma1 = ( X1 - ones(n1,1)*mu1' )' *  ( X1 - ones(n1,1)*mu1' ) / n1;
        Sigma0 = ( X0 - ones(n0,1)*mu0' )' *  ( X0 - ones(n0,1)*mu0' ) / n0;
        
        Sigma = diag(diag(n0/n * Sigma0 + n1/n * Sigma1));
        
        % computing discriminative affine function
        w_LDA  = Sigma \ ( mu1 - mu0 );
        b_LDA = - .5 * ( mu1 + mu0 )' * w_LDA + log ( n1 / n0 );
        
        
        % compute error rates
        error_train_LDA_diag(irep,it) = mean( abs( ( sign( X * w_LDA + b_LDA ) + 1 ) / 2 - y ) )
        error_test_LDA_diag(irep,it) = mean( abs( ( sign( X_test * w_LDA + b_LDA ) + 1 ) / 2 - y_test ) )
        
        
        % LOGISTIC REGRESSION
        Xint = [ X ones(n,1) ];   % add constant variable to x
        theta = zeros(d+1,1);
        normgradient = Inf;
        k=1;
        while (normgradient > 1e-10 & k < 100)
            eta = Xint * theta;
            mu = 1./(1+exp(-eta));
            gradient = Xint'*(y-mu) - 1e-12 * theta;
            H = - Xint' * ( repmat( mu.*(1-mu) , 1 , d+1 ) .* Xint ) - 1e-12 * eye(d+1);
            theta = theta - H \ gradient;
            if any(isnan(theta)), pause; end
            normgradient = norm(gradient);
            fprintf('k = %d - norm of gradient = %e\n',k,normgradient);
            k=k+1;
        end
        theta_LOGISTIC = theta;
        w_LOGISTIC = theta(1:d);
        b_LOGISTIC = theta(d+1);
        % compute error rates
        error_train_LOGISTIC(irep,it) = mean( abs( ( sign( X * w_LOGISTIC + b_LOGISTIC ) + 1 ) / 2 - y ) )
        error_test_LOGISTIC(irep,it) = mean( abs( ( sign( X_test * w_LOGISTIC + b_LOGISTIC ) + 1 ) / 2 - y_test ) )
        
    end
end
subplot(1,2,2);
plot(ts,mean(error_test_LDA_diag,1),'r','linewidth',2); hold on;
plot(ts,mean(error_test_LDA,1),'b','linewidth',2); hold on;
plot(ts,mean(error_test_LOGISTIC,1),'k','linewidth',2); hold off;
legend('LDA Diag','LDA','logistic','location','southeast');
set(gca,'fontsize',16)
xlabel('lack of independence');
ylabel('error rate')
axis([0 1 0 .23]);
title('d = 50, n = 500', 'fontweight','normal');

try
    print('-depsc', 'compare_LDA_Logistic.eps');
    close(ccc)
catch
    disp('missing figure file')
end

