%%

k = [1,3, 5, 7, 9];
c = [1, 3, 5, 7, 9];
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

R = wblrnd(k(5),c(3),[1000,1]);
m = mean(R)
v = var(R)
kVar = sqrt(v)/m;

Fx = 1 - exp((-x/c(3)).^k(5))
 = x.*
