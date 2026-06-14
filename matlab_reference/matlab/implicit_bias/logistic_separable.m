clear all
seed=1;
randn('state',seed);
rand('state',seed);

try
    ccc=openfig('logistic_separable.fig');
catch
    disp('missing figure file')
end




% separable data
d = 2;
n = 100;
y = [ones(n/2,1); -ones(n/2,1) ];
X = [-ones(n/2,1) -ones(n/2,1) ; ones(n/2,1) ones(n/2,1) ] + randn(n,2)/2;
plot(X(:,1),X(:,2))


X = [X, ones(n,1) * 1];
d = d + 1;
L = max(eig(X'*X)/n);
w = randn(d,1);
w(end) = 0;


% GD
maxiter  = 100000;
for iter=1:maxiter
    ws(:,iter) = w;
    functionerror(iter) = mean(    y .* (  X*w )  < 0 ) ;
    functionval(iter) = mean(  log( 1 + exp(- y .* (  X*w ) ) ) ) ;
    temp = ( X * w ) .* y;
    grad =   1/n * X' *  ( y .* ( - 1./( 1 + exp(temp) ) ) ) ;
    w = w - 1/(L) * grad;
end

itertoplot = round( 2.^(.5:.5:9))
for i=1:length(itertoplot)
    it=itertoplot(i);
    subplot(3,6,i);
    plot(X(find(y==1),1),X(find(y==1),2),'bx','markersize',10); hold on;
    plot(X(find(y==-1),1),X(find(y==-1),2),'rx','markersize',10);
    Xright = 10;
    Xleft = -10;
    
    Yleft = ( ws(3,it) - Xleft*ws(1,it) ) / ws(2,it);
    Yright = ( ws(3,it) - Xright*ws(1,it) ) / ws(2,it);
    plot([Xleft, Xright], [Yleft, Yright],'k','linewidth',2);
    
    Yleft = ( ws(3,end) - Xleft*ws(1,end) ) / ws(2,end);
    Yright = ( ws(3,end) - Xright*ws(1,end) ) / ws(2,end);
    plot([Xleft, Xright], [Yleft, Yright],'k:','linewidth',2);
    
    hold off
    axis([-4 4 -4 4]);
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])
set(gca,'fontsize',20);
    title(sprintf('t = %d',it),'FontWeight','normal');
end

try
    print('-depsc', 'logistic_separable.eps');
    close(ccc)
catch
    disp('missing figure file')
end

