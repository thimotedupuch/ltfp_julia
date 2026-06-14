clear all

try
    ccc=openfig('mmse.fig');
catch
    disp('missing figure file')
end



% spike and slab



y = -5:.015:5;


subplot(1,3,1);
alpha = .95;
tau = sqrt(1/(1-alpha)); % make it unit variance


% noise variances
sigmas = [ 1/8 1/2 1 2];

for isigma = 1:length(sigmas)
    sigma = sigmas(isigma);   
    py = alpha / sigma /sqrt(2*pi) * exp( - y.^2 / 2/sigma^2) +  (1-alpha) / sqrt(sigma^2+tau^2) /sqrt(2*pi) * exp( - y.^2 / 2/(sigma^2+tau^2));    
    MMSE(:,isigma) = y(2:end) + sigma^2 * ( log(py(2:end)) - log(py(1:end-1)) )/(y(2)-y(1));
end
plot(y(2:end),MMSE,'linewidth',2)
set(gca,'fontsize',18)
legend('1/8','1/2','1','2','Location','NorthWest');
title('Spike and Slab - MMSE','FontWeight','normal');
xlabel('y');
ylabel('MMSE(y)');





subplot(1,3,2);
lambda = sqrt(2); % make it unit variance
for isigma = 1:length(sigmas)
    sigma = sigmas(isigma);    
    py = lambda/4 * exp( lambda^2 * sigma^2 / 2 - lambda * y) .* ( 1 - erf((lambda *sigma -y/sigma)/sqrt(2))) +   lambda/4 * exp( lambda^2 * sigma^2 / 2 + lambda *  y) .* ( 1 - erf((lambda*sigma+y/sigma)/sqrt(2)))    
    MMSE(:,isigma) = y(2:end) + sigma^2 * ( log(py(2:end)) - log(py(1:end-1)) )/(y(2)-y(1));
    MAP(:,isigma) = sign(y(2:end)) .* max( abs(y(2:end)) - lambda * sigma^2 , 0);   
end
plot(y(2:end),MMSE,'linewidth',2); 
set(gca,'fontsize',18)
legend('1/8','1/2','1','2','Location','NorthWest');
title('Laplace - MMSE','FontWeight','normal');
xlabel('y');
ylabel('MMSE(y)');


subplot(1,3,3);
plot(y(2:end),MAP,'linewidth',2);  
set(gca,'fontsize',18)
legend('1/8','1/2','1','2','Location','NorthWest');
title('Laplace - MAP','FontWeight','normal');
xlabel('y');
ylabel('MAP(y)');


try
    print('-depsc', 'MMSE.eps');
    close(ccc)
catch
    disp('missing figure file')
end




% alpha = 0;
% tau = sqrt(1/(1-alpha));
% 
% 
% for isigma = 1:length(sigmas)
%     sigma = sigmas(isigma);
%     
%     py = alpha / sigma /sqrt(2*pi) * exp( - y.^2 / 2/sigma^2) +  (1-alpha) / sqrt(sigma^2+tau^2) /sqrt(2*pi) * exp( - y.^2 / 2/(sigma^2+tau^2));
%     
%     
%     MMSE(:,isigma) = y(2:end) + sigma^2 * ( log(py(2:end)) - log(py(1:end-1)) )/(y(2)-y(1));
% end
% plot(y(2:end),MMSE,'linewidth',2)
% set(gca,'fontsize',16)
% legend('1/8','1/2','1','2','Location','NorthWest');
% title('Gaussian','FontWeight','normal');
% xlabel('y');
% ylabel('MMSE(y)');
% 
% 
% 
