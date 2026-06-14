clear all
seed = 1;
rand('state',seed);
randn('state',seed);
addpath ..

try
ccc=openfig('robust_regression.fig');
catch
    disp('missing figure file')
end


n = 100;
x = rand(n,1);
xtest=(0:.001:1)';
std_noise = .1;
y = sin(x.^2 * 6*pi) + std_noise * randn(n,1).^5;
ytest = sin(xtest.^2 * 6*pi);
plot(x,y,'x'); hold on
plot(xtest,ytest,'r-'); hold off

alpha_kernel = 10;
K = exp(-alpha_kernel*sq_dist(x',x'));
Ktest = exp(-alpha_kernel*sq_dist(xtest',x'));
lambdas = 10.^[-12:.5:1];

% square loss
for ilambda = 1:length(lambdas)
    
    
    lambda = lambdas(ilambda);
    
    
    alpha = ( K + n * lambda  *eye(n)) \ y;
    ypredtest = Ktest * alpha;
    
    
    square_error(ilambda) = mean( (ypredtest - ytest).^2)
    robust_error(ilambda) = mean(  1- exp( -(ypredtest - ytest).^2) )
    
end

[a,ilambda] = min(robust_error);
lambda = lambdas(ilambda);
alpha = ( K + n * lambda  *eye(n)) \ y;
ypredtest = Ktest * alpha;

subplot(1,2,1)
plot(xtest,ytest,'r-','linewidth',2); hold on
plot(xtest,ypredtest,'k-','linewidth',2);
plot(x,y,'bx','markersize',10); hold off
xlabel('x');
ylabel('y');
set(gca,'fontsize',16)
legend('target','prediction')
title('square loss', 'fontweight','normal');


% robust loss
% square loss
for ilambda = 1:length(lambdas)
    lambda = lambdas(ilambda);
    
    
    for itest =1:length(ytest)
        
        temp =   ( K + n * lambda  *eye(n)) \  Ktest(itest,:)';
        z = -5:.01:5;
        
        [a,b] = min( ( 1-exp(-(z-y).^2) )' *  temp );
        
        ypredtest(itest) =z(b);
    end
    square_error_rob(ilambda) = mean( (ypredtest - ytest).^2)
    robust_error_rob(ilambda) = mean(  1- exp( -(ypredtest - ytest).^2) )
    
end

[a,ilambda] = min(robust_error);
lambda = lambdas(ilambda);


for itest =1:length(ytest)
    
    temp =   ( K + n * lambda  *eye(n)) \  Ktest(itest,:)';
    z = -5:.005:5;
    
    [a,b] = min( ( 1-exp(-(z-y).^2) )' *  temp );
    
    ypredtest(itest) =z(b);
end

subplot(1,2,2)
plot(xtest,ytest,'r-','linewidth',2); hold on
plot(xtest,ypredtest,'k-','linewidth',2);
plot(x,y,'bx','markersize',10); hold off
xlabel('x');
ylabel('y');
set(gca,'fontsize',16)
legend('target','prediction')
title('robust loss', 'fontweight','normal');

try
print('-depsc', 'robust_regression.eps');
close(ccc)
catch
    disp('missing figure file')
end

