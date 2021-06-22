function ADMM_graph_result(hist, k, MLSFji)
figure(19)
pp1 = semilogy(1:(k-1),[hist.primal],'LineWidth',2,'color','k');
hold on
pp2 = semilogy(1:(k-1),[hist.dual],'--','LineWidth',2,'color','k');
hold off
grid on
xlabel('iterations','FontSize',12)
ylabel('residuals','FontSize',12)
set(gca,'FontSize',16)
leg0 = legend([{"primal residual"}, {"dual residual"}]);

es = [];
eb = [];
eaverage = [];
zglob = [];
sellprice = [];
buyprice = [];

f20 = figure(20);
% f20.Position = [0, 0, 650, 650];
NumIterations = 150;
yaxislimit = 100;

for i=1:k-1
    es(i) = hist(i).Es(4,2);
    eb(i) = hist(i).Eb(4,2);
    eaverage(i) = hist(i).Eaverage(4,2);
    zglob(i) = hist(i).zglob(4,2);
    sellprice(i) = hist(i).Lambda_s(4,2);
    buyprice(i) = hist(i).Lambda_b(4,2);
end
subplot(2,2,1)
pp1 = plot(1:(k-1),es,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),eb,'-','LineWidth',1.5, 'color',[1 0.5 0]);
pp3 = plot(1:(k-1),eaverage,'-','LineWidth',1.5);
pp4 = plot(1:(k-1),zglob,'--','LineWidth',1.5,'color','k');
hold off
grid on
leg1 = legend([{"$e^s_{24}$"}, {"$e^b_{24}$"}, {"$\bar{e}_{24}$"}, {"$z_{24}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none','numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
xticks([0 50 100 150])
yticks([0 25 50 75 100])
ylabel('Energy(kWh)','FontSize',12)
set(gca,'FontSize',16)
title('seller2-buyer4')

subplot(2,2,2)
pp1 = plot(1:(k-1),sellprice,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',1.5, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{24}$"}, {"$\lambda^b_{24}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([2 8])
xticks([0 50 100 150])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller2-buyer4')

for i=1:k-1
    es(i) = hist(i).Es(1,5);
    eb(i) = hist(i).Eb(1,5);
    eaverage(i) = hist(i).Eaverage(1,5);
    zglob(i) = hist(i).zglob(1,5);
    sellprice(i) = hist(i).Lambda_s(1,5);
    buyprice(i) = hist(i).Lambda_b(1,5);
end
subplot(2,2,3)
pp1 = plot(1:(k-1),es,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),eb,'-','LineWidth',1.5, 'color',[1 0.5 0]);
pp3 = plot(1:(k-1),eaverage,'-','LineWidth',1.5);
pp4 = plot(1:(k-1),zglob,'--','LineWidth',1.5,'color','k');
hold off
grid on
leg1 = legend([{"$e^s_{51}$"}, {"$e^b_{51}$"}, {"$\bar{e}_{51}$"}, {"$z_{51}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none','numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
yticks([0 25 50 75 100])
xticks([0 50 100 150])
ylabel('Energy(kWh)','FontSize',12)
set(gca,'FontSize',16)
title('seller5-buyer1')

subplot(2,2,4)
pp1 = plot(1:(k-1),sellprice,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',1.5, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{51}$"}, {"$\lambda^b_{51}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([2 8])
xticks([0 50 100 150])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller5-buyer1')

f21 = figure(21);
% f21.Position = [0, 0, 650, 650];
NumIterations = 75;
yaxislimit = 100;

for i=1:k-1
    es(i) = hist(i).Es(3,2);
    eb(i) = hist(i).Eb(3,2);
    eaverage(i) = hist(i).Eaverage(3,2);
    zglob(i) = hist(i).zglob(3,2);
    sellprice(i) = hist(i).Lambda_s(3,2);
    buyprice(i) = hist(i).Lambda_b(3,2);
end
subplot(2,2,1)
pp1 = plot(1:(k-1),es,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),eb,'-','LineWidth',1.5, 'color',[1 0.5 0]);
pp3 = plot(1:(k-1),eaverage,'-','LineWidth',1.5);
pp4 = plot(1:(k-1),zglob,'--','LineWidth',1.5,'color','k');
hold off
grid on
leg1 = legend([{"$e^s_{21}$"}, {"$e^b_{21}$"}, {"$\bar{e}_{21}$"}, {"$z_{21}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none','numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
xticks([0 25 50 75])
yticks([0 25 50 75 100])
ylabel('Energy(kWh)','FontSize',12)
set(gca,'FontSize',16)
title('seller2-buyer1')

subplot(2,2,2)
pp1 = plot(1:(k-1),sellprice,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',1.5, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{21}$"}, {"$\lambda^b_{21}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([2 8])
xticks([0 25 50 75])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller2-buyer1')

for i=1:k-1
    es(i) = hist(i).Es(2,6);
    eb(i) = hist(i).Eb(2,6);
    eaverage(i) = hist(i).Eaverage(2,6);
    zglob(i) = hist(i).zglob(2,6);
    sellprice(i) = hist(i).Lambda_s(2,6);
    buyprice(i) = hist(i).Lambda_b(2,6);
end
subplot(2,2,3)
pp1 = plot(1:(k-1),es,'LineWidth',2, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),eb,'-','LineWidth',2, 'color',[1 0.5 0]);
pp3 = plot(1:(k-1),eaverage,'-','LineWidth',2);
pp4 = plot(1:(k-1),zglob,'--','LineWidth',2,'color','k');
hold off
grid on
leg1 = legend([{"$e^s_{15}$"}, {"$e^b_{15}$"}, {"$\bar{e}_{15}$"}, {"$z_{15}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none','numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
yticks([0 25 50 75 100])
xticks([0 25 50 75])
ylabel('Energy(kWh)','FontSize',12)
set(gca,'FontSize',16)
title('seller1-buyer5')

subplot(2,2,4)
pp1 = plot(1:(k-1),sellprice,'LineWidth',2, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',2, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{15}$"}, {"$\lambda^b_{15}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([2 8])
xticks([0 25 50 75])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller1-buyer5')

f22 = figure(22);
% f21.Position = [0, 0, 650, 650];
NumIterations = 100;
yaxislimit = 100;

for i=1:k-1
    zglob(i) = hist(i).zglob(3,5);
    loss(i) = zglob(i)*(MLSFji(3,5));
    sellprice(i) = hist(i).Lambda_s(3,5);
    buyprice(i) = hist(i).Lambda_b(3,5);
end
subplot(2,2,1)
yyaxis right
pp2 = area(1:(k-1),loss,'FaceColor',[1,0,0],'FaceAlpha',0.4,'EdgeColor','None');
ylim([-4 4])
xlim([0 NumIterations])
yyaxis left
pp1 = plot(1:(k-1),zglob,'--','LineWidth',2,'color','k');
ylim([0 yaxislimit])
xlim([0 NumIterations])
xticks([0 25 50 75 100])
yticks([0 25 50 75 100])
ylabel('Energy(kWh)','FontSize',12)
grid on
leg1 = legend([{"$z_{24}$"},{"$\tau_{24}z_{24}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none');
set(gca,'FontSize',16)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
title('seller2-buyer4')

subplot(2,2,2)
pp1 = plot(1:(k-1),sellprice,'LineWidth',2, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',2, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{24}$"}, {"$\lambda^b_{24}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([4 5])
xticks([0 25 50 75 100])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller2-buyer4')

for i=1:k-1
    zglob(i) = hist(i).zglob(4,4);
    loss(i) = zglob(i)*(MLSFji(4,4));
    sellprice(i) = hist(i).Lambda_s(4,4);
    buyprice(i) = hist(i).Lambda_b(4,4);
end
subplot(2,2,3)
yyaxis right
pp2 = area(1:(k-1),loss,'FaceColor',[1,0,0],'FaceAlpha',0.4,'EdgeColor','None');
ylim([-2 2])
xlim([0 NumIterations])
yyaxis left
pp1 = plot(1:(k-1),zglob,'--','LineWidth',2,'color','k');
ylim([0 yaxislimit])
xlim([0 NumIterations])
xticks([0 25 50 75 100])
yticks([0 25 50 75 100])
ylabel('Energy(kWh)','FontSize',12)
grid on
leg1 = legend([{"$z_{33}$"},{"$\tau_{33}z_{33}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'color','none');
set(gca,'FontSize',16)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
title('seller3-buyer3')

subplot(2,2,4)
pp1 = plot(1:(k-1),sellprice,'LineWidth',2, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'--','LineWidth',2, 'color',[1 0.5 0]);
hold off
grid on
leg2 = legend([{"$\lambda^s_{33}$"}, {"$\lambda^b_{33}$"}]);
set(leg2,'Interpreter','latex','FontSize',16,'color','none');
xlim([0 NumIterations])
ylim([4 5])
xticks([0 25 50 75 100])
ylabel("Price("+char(0162)+"/kWh)",'FontSize',12)
set(gca,'FontSize',16)
title('seller3-buyer3')
end