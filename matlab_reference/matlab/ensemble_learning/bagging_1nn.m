clear all
seed = 1;
randn('state',seed);
rand('state',seed);

xtest = (0:.001:1)';
fxtest = sin(2*pi*xtest);

std_noise = .35;
n = 200;
x = ( 0:(n-1) )'  / n;
x = sort(x);
ybeforenoise = sin(2*pi*x);
y = sin(2*pi*x) + std_noise * randn(n,1);

ntest = length(xtest);

ms = unique(round(2.^[0:1:7]));
nrep = 100;

try
    ccc=openfig('bagging_1nn.fig');
catch
    disp('missing figure file')
end



for im = 1:length(ms)
    m = ms(im);
    nloc = round(n / m);
    
    for irep=1:nrep
        ind = randperm(n);
        ind = ind(1:nloc);
        dists = sq_dist(x(ind)',xtest');
        [a,b] = min(dists,[],1);
        ytestpred(irep,:) = y(ind(b));
    end
    
    subplot(2,4,im);
    plot(xtest,mean(ytestpred,1),'b','linewidth',2); hold on;
    plot(xtest,fxtest,'r','linewidth',2); hold on;
    plot(x,y,'kx'); hold off;
    title(sprintf('n/s = %d - error = %1.2e',m,mean( (mean(ytestpred,1)-fxtest').^2)),'fontweight','normal');
    set(gca,'fontsize',18);
    legend('bagging','target');
    axis([0 1 -1.5 1.6]);
end
try
    print('-depsc', 'bagging_1nn.eps');
    close(ccc)
catch
    disp('missing figure file')
end


