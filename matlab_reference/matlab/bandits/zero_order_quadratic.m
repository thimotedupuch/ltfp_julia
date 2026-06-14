clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('zeroorder.fig');
catch
    disp('missing figure file')
end


% fixed matrix with planted singular values
d = 10;
H = randn(d,d);
[u,s,v] = svd(H);
H = u * diag(1./(1:d)) * u';
% H = u * diag(1./(1:d).^2) * u';
Hsqrt = sqrtm(H);
L = max(eig(H));

theta_ast = randn(d,1);
sigma = .01;
maxiter = 100000;
deltas = [ 0.01 0.1 10]
nrep = 10

val_zero = zeros(length(deltas),maxiter,nrep);
val_zero_ave = zeros(length(deltas),maxiter,nrep);
for idelta = 1:length(deltas);
    delta = deltas(idelta);
    
    for irep=1:nrep
        w = zeros(d,1);
        wave = zeros(d,1);
        for iter=1:maxiter
            val_zero(idelta,iter,irep) = .5 * (w - theta_ast)'*H*(w-theta_ast);
            val_zero_ave(idelta,iter,irep) = .5 * (wave - theta_ast)'*H*(wave-theta_ast);
            
            
            
            z = randn(d,1);
            temp1 = .5 * (w + delta * z - theta_ast)'*H*(w + delta * z-theta_ast);
            temp2 = .5 * (w - theta_ast)'*H*(w-theta_ast);
            
            w = w - 1/4/(L)/d * ( temp1 - temp2 +randn(1) * sigma) / delta * z ;
            wave = ( 1 - 1/iter) * wave + w/iter;
        end
    end
end
subplot(2,3,1);
plot( log10(mean(val_zero(1,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(1,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(1),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging');

subplot(2,3,2);
plot( log10(mean(val_zero(2,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(2,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(2),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging');

subplot(2,3,3);
plot( log10(mean(val_zero(3,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(3,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(3),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging','location','southeast');





sigma = .1;
maxiter = 100000;
deltas = [ 0.01 0.1 10]
nrep = 10

val_zero = zeros(length(deltas),maxiter,nrep);
val_zero_ave = zeros(length(deltas),maxiter,nrep);
for idelta = 1:length(deltas);
    delta = deltas(idelta);
    
    for irep=1:nrep
        w = zeros(d,1);
        wave = zeros(d,1);
        for iter=1:maxiter
            val_zero(idelta,iter,irep) = .5 * (w - theta_ast)'*H*(w-theta_ast);
            val_zero_ave(idelta,iter,irep) = .5 * (wave - theta_ast)'*H*(wave-theta_ast);
            
            
            
            z = randn(d,1);
            temp1 = .5 * (w + delta * z - theta_ast)'*H*(w + delta * z-theta_ast);
            temp2 = .5 * (w - theta_ast)'*H*(w-theta_ast);
            
            w = w - 1/4/(L)/d * ( temp1 - temp2 +randn(1) * sigma) / delta * z ;
            wave = ( 1 - 1/iter) * wave + w/iter;
        end
    end
end
subplot(2,3,4);
plot( log10(mean(val_zero(1,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(1,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(1),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging','location','southeast');

subplot(2,3,5);
plot( log10(mean(val_zero(2,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(2,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(2),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging');

subplot(2,3,6);
plot( log10(mean(val_zero(3,:,:),3)),'b','linewidth',2); hold on;
plot( log10(mean(val_zero_ave(3,:,:),3)),'r','linewidth',2); hold off
axis([0 maxiter -5 1.5]);
title(sprintf('\\delta = %1.2f - \\sigma = %1.2f',deltas(3),sigma),'fontweight','normal')
set(gca,'fontsize',18)
xlabel('# iterations');
ylabel('F(\theta)-F(\theta_*)');
legend('no averaging','averaging','location','southeast');

try
    print('-depsc', 'zeroorder.eps');
    close(ccc)
catch
    disp('missing figure file')
end

