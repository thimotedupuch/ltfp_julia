clear all
seed= 4 ;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('ucb.fig');
catch
    disp('missing figure file')
end

nrep = 1;
k = 10;
maxt = 10000;
means = rand(k,1);
rewards = zeros(maxt,nrep);
for irep=1:nrep
    % initializing means
    rewards_per_arm = double(rand(k,1) < means);
    rewards(1:k,irep) = rewards_per_arm;
    ns = ones(k,1);
    rho = 2*log(k);
    upper_bounds(:,k) = rewards_per_arm ./ ns + sqrt(.5*rho./ns);
    lower_bounds(:,k) = rewards_per_arm ./ ns - sqrt(.5*rho./ns);
    
    for t=k+1:maxt
        [a,arm] = max( upper_bounds(:,t-1) );
        reward = double( rand < means(arm));
        rewards(t,irep) = reward;
        rewards_per_arm(arm) = rewards_per_arm(arm)+reward;
        ns(arm) = ns(arm)+1;
        rho = 2*log(t);
    upper_bounds(:,t) = rewards_per_arm ./ ns + sqrt(.5*rho./ns);
    lower_bounds(:,t) = rewards_per_arm ./ ns - sqrt(.5*rho./ns);
        
    end
end
subplot(1,2,1);
plot(0,0,'r','linewidth',2); hold on; plot(0,0,'b','linewidth',2);  

plot(lower_bounds','.r-','linewidth',2); hold on; plot(upper_bounds','.b-','linewidth',2); hold off
set(gca,'fontsize',18)
xlabel('t');
legend('lower','upper');
ylabel('mean values');
title('confidence bounds','fontweight','normal');




nrep = 100;
k = 10;
maxt = 10000;
means = rand(k,1);
rewards = zeros(maxt,nrep);
for irep=1:nrep
    % initializing means
    rewards_per_arm = double(rand(k,1) < means);
    rewards(1:k,irep) = rewards_per_arm;
    ns = ones(k,1);
    rho = 2*log(k);
    upper_bounds  = rewards_per_arm ./ ns + sqrt(.5*rho./ns);
     
    for t=k+1:maxt
        [a,arm] = max( upper_bounds  );
        reward = double( rand < means(arm));
        rewards(t,irep) = reward;
        rewards_per_arm(arm) = rewards_per_arm(arm)+reward;
        ns(arm) = ns(arm)+1;
        rho = 2*log(t);
    upper_bounds = rewards_per_arm ./ ns + sqrt(.5*rho./ns);
         
    end
end

  
 perf = mean(cumsum( max(means) - rewards(:,:)) ,2);
 stdperf = std(cumsum( max(means) - rewards(:,:)),[],2);
 
 [a,b] = affine_fit(log(maxt/10:maxt), perf(maxt/10:maxt))
 
subplot(1,2,2);
plot(1:maxt,perf,'r','linewidth',2);
hold on
plot(1:maxt, log(1:maxt)*a+b,'b','linewidth',2);
% plot(1:maxt,perf+stdperf,'r:','linewidth',2);
% plot(1:maxt,perf-stdperf,'r:','linewidth',2);
hold off

legend('observation','fit: a+b*log(t)','location','southeast')
set(gca,'fontsize',18)
xlabel('t');
ylabel('regret');
axis([1 maxt 0 220])
 title('regret','fontweight','normal');
try
    print('-depsc', 'ucb_matlab.eps');
    close(ccc)
catch
    disp('missing figure file')
end

