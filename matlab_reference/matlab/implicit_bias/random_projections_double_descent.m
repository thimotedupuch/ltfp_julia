clear all
seed= 0 ;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('doubledescent_top.fig');
catch
    disp('missing figure file')
end



n = 200;
d = 400;
temp = randn(d,d); temp = temp + temp';
[u,e] = eig(temp);
Sigma = u*diag(1./( (1:d) ))*u';
Sigma = Sigma / trace(Sigma);

Sigmasqrt = sqrtm(Sigma);
nrep = 20;
ms = 0:2:2*d;
ms( find(abs(ms-n)<5) ) = []; % remove the m's which are too close to n

empRvar = Inf * ones(length(ms),nrep);
empRbias = Inf * ones(length(ms),nrep);


% single prediction problem
theta_ast =   randn(d,1)  ;
theta_ast = theta_ast / sqrt( theta_ast'*Sigma*theta_ast);



for irep=1:nrep
    
    % X is averaged out
    
    X = sign(randn(n,d)) * Sigmasqrt;
    invK = inv(X*X');
    Pi = X' * invK * X;
    
    
    %  S and epsilon are averaged out
    
    epsilon = sign(randn(n,1));
    Sfull = sign(randn(d,max(ms)));
    
    
    for im = 1:length(ms);
        m = ms(im);
        S = Sfull(:,1:m);
        
        if m < n
            
            invSXXS = inv( S'*X'*X*S );
            thetabias =  S * invSXXS * S' * X' * X * theta_ast;
            thetavar =  S * invSXXS * S' * X' * epsilon;
            
            empRvar(im,irep) = (thetavar)'*Sigma*(thetavar);
            empRbias(im,irep) = (thetabias-theta_ast)'*Sigma*(thetabias-theta_ast);
            
        elseif  m > n
            
            
            invXSSX = inv(X*S*S'*X');
            thetabias =   S*S'*X'* invXSSX * X * theta_ast;
            thetavar =  S*S'*X'* invXSSX * epsilon;
            
            empRvar(im,irep) = (thetavar)'*Sigma*(thetavar);
            empRbias(im,irep) = (thetabias-theta_ast)'*Sigma*(thetabias-theta_ast);
            
        else
            invXSSX = inv(X*S*S'*X');
            thetabias =   S*S'*X'* invXSSX * X * theta_ast;
            thetavar =  S*S'*X'* invXSSX * epsilon;
            
            empRvar(im,irep) = (thetavar)'*Sigma*(thetavar);
            empRbias(im,irep) = (thetabias-theta_ast)'*Sigma*(thetabias-theta_ast);
            
        end
        
    end
end

std_noise = 1/2;

plot(ms,mean(empRbias,2)+ std_noise^2 * mean(empRvar,2),'b','linewidth',2); hold on
plot(ms,mean(empRbias,2)+ std_noise^2 * mean(empRvar,2) + std(empRbias +  std_noise^2 * empRvar, [],2),'b:','linewidth',2);
temp = mean(empRbias,2)+ std_noise^2 * mean(empRvar,2) - std(empRbias +  std_noise^2 * empRvar,[],2);
temp = ( temp < 0 ) * 1e12 + temp; % removes the really unstable ones that we don't see anyway
plot(ms,temp,'b:','linewidth',2); hold off
axis([ 0 max(ms) 0 4])
set(gca,'fontsize',24);
xlabel('m');
ylabel('excess risks');

try
    print('-depsc', 'doubledescent_randomfeature.eps');
    close(ccc)
catch
    disp('missing figure file')
end

