%%

k = [1,3, 5, 7, 9];
c = [1, 3, 5, 7, 9];
n = 1000;
B = 500;
close all;
for i=1:5
    for j=1:5
        x = linspace(0,10,1000);
        fx = k(i)/c(j)*(x/c(j)).^(k(i)-1).*exp((-x/c(j)).^k(i));
        % figure(i)
        % %subplot(5,1,i);
        % plot(x,fx); hold on;
    end
end

R1 = wblrnd(k(5),c(3),[n,1]);
    m1 = mean(R1)
    v1 = var(R1)
    kVar = sqrt(v1)/m1
    
R = zeros(n,B);
ERT = zeros(1,B);
kMM =  zeros(1,B);
cMM = zeros(1,B);

tic
for a=1:500
    R(:,a) = wblrnd(k(5),c(3),[n,1]);
    m = mean(R(:,a));
    v = var(R(:,a));
    kVar = sqrt(v)/m;
    
    
    m2 = (R(:,a)'*R(:,a))/n;
    m1 = sum(R(:,a))/n;
    m1Square  = m1^2;
    valkIdeal = m2/m1Square;
    %Fx = 1 - exp((-x/c(3)).^k(5))
    % = x.*
    
    kI = 10;
    err = valkIdeal + 1;
    i = 1;
    while i < 5000 && abs(err)>1e-5 % petite acceleration de Steffensen AIIIGGGHT :D
        
        FkI =  2*kI*gamma(2/kI)/(gamma(1/kI))^2 - valkIdeal;%gamma(2/kI+1)/(gamma(1/kI+1))^2 - valkIdeal;
        fkI = kI - FkI/4;
        FfKI = 2*fkI*gamma(2/fkI)/(gamma(1/fkI))^2 - valkIdeal;%gamma(2/fkI+1)/(gamma(1/fkI+1))^2 - valkIdeal;
        ff = fkI - FfKI/2;
        kII = kI - (fkI-kI)^2/(ff-2*fkI+kI);
        err = kII - kI;
        kI = kII;
        i = i +1 ;
        
    end
    
    kMM(1,a) = kI;
    cMM(1,a) = kMM(1,a)*m1/gamma(1/kMM(1,a));
    
    ERT(1,a) = (kMM(1,a) - k(5))^2 + (cMM(1,a)-c(3))^2;
    
end
toc
ERTmean = mean(ERT)
kMMmean = mean(kMM)
cMMmean = mean(cMM)

ERTvar = var(ERT)
kMMvar = var(kMM)
cMMvar = var(cMM)

figure(1)
boxplot(ERT)
figure(2)
boxplot(kMM)
figure(3)
boxplot(cMM)

figure(4)
hist(ERT)
figure(5)
hist(kMM)
figure(6)
hist(cMM)

