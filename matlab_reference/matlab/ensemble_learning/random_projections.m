clear all
seed = 1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('random_projection.fig');
catch
    disp('missing figure file')
end



n = 100;
ntest = 100;
dx = 20;
nrep = 100;
d = dx*(dx+1)/2+dx+1;

for irep=1:nrep
    
    X = randn(n,dx);
    Xtest = X; %randn(ntest,dx);
    Phi = zeros(n,d);
    Phi(:,1) = 1;
    Phi(:,2:dx+1) = X;
    k = dx+2;
    for i=1:dx
        for j=1:i
            Phi(:,k)   = X(:,i) .* X(:,j);
            k = k + 1;
        end
    end
    
    Phitest = zeros(ntest,d);
    Phitest(:,1) = 1;
    Phitest(:,2:dx+1) = Xtest;
    k = dx+2;
    for i=1:dx
        for j=1:i
            Phitest(:,k)   = Xtest(:,i) .* Xtest(:,j);
            k = k + 1;
        end
    end
    
    
    std_noise = 2;
    y = X(:,1).^2 + std_noise * randn(n,1);
    ytest =  Xtest(:,1).^2 + std_noise * randn(ntest,1); ;
    
    % ridge regression
    lambdas = 10.^[4:-.5:-8];
    for ilambda = 1:length(lambdas)
        lambda = lambdas(ilambda);
        
        theta = ( Phi'*Phi + n*lambda*eye(d) ) \ ( Phi'* y);
        ypred = Phi * theta;
        ypredtest = Phitest * theta;
        
        train_error_ridge(ilambda,irep) = mean( (ypred-y).^2 );
        test_error_ridge(ilambda,irep) = mean( (ypredtest-ytest).^2 );
    end
    
    
    % random projection
    ss = 0:5:100;
    ms = [1 10 100];
    for imm = 1:length(ms)
        m = ms(imm);
        
        for im=1:m
            Ss{im} = randn(d,max(ss));
        end
        for is = 1:length(ss)
            s = ss(is);
            thetas=[];
            for im=1:m
                S = Ss{im}(:,1:s);
                thetas(:,im) = S * ( ( S'*Phi'*Phi*S) \ ( S'*Phi'* y) );
            end
            theta = mean(thetas,2);
            ypred = Phi*theta;
            ypredtest = Phitest*theta;
            train_error_rp(is,imm,irep) = mean( (ypred-y).^2 );
            test_error_rp(is,imm,irep) = mean( (ypredtest-ytest).^2 );
        end
    end
end

% plotting
subplot(2,2,1)
errorbar(-log(lambdas),mean(train_error_ridge,2),std(train_error_ridge,[],2), std(train_error_ridge,[],2),'b','linewidth',2); hold on;
errorbar(-log(lambdas),mean(test_error_ridge,2), std(test_error_ridge,[],2), std(test_error_ridge,[],2),'r','linewidth',2); hold off;
legend('train','test');
set(gca,'fontsize',20)
title('Ridge regression','fontweight','normal');
xlabel('-log(\lambda)')
ylabel('errors')
axis([min(-log(lambdas)) max(-log(lambdas)) 0 10])


for im=1:3
subplot(2,2,1+im)
errorbar(ss,mean(train_error_rp(:,im,:),3),std(train_error_rp(:,im,:),[],3), std(train_error_rp(:,im,:),[],3),'b','linewidth',2); hold on;
errorbar(ss,mean(test_error_rp(:,im,:),3), std(test_error_rp(:,im,:),[],3), std(test_error_rp(:,im,:),[],3),'r','linewidth',2); 
hold off
legend('train','test');
set(gca,'fontsize',20)
axis([min(ss) max(ss) 0 10])
title(sprintf('Random projections - m=%d',ms(im)),'fontweight','normal');
xlabel('s')
ylabel('errors')
end



try
    print('-depsc', 'random_projection.eps');
    close(ccc)
catch
    disp('missing figure file')
end

