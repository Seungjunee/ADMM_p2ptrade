function ADMM_graph_result(hist, k)
figure(19)
pp1 = semilogy(1:(k-1),[hist.primal],'LineWidth',2,'color','k');
hold on
pp2 = semilogy(1:(k-1),[hist.dual],'--','LineWidth',2,'color','k');
hold off
grid on
xlabel('iterations','FontSize',14)
ylabel('residuals','FontSize',14)
set(gca,'FontSize',14)
leg0 = legend([{"primal residual"}, {"dual residual"}]);

es = [];
eb = [];
eaverage = [];
zglob = [];
sellprice = [];
buyprice = [];

f20 = figure(20);
% set(f20,'Position',[0 0 600 400])
% f20.Position = [0, 0, 650, 650];
NumIterations = 150;
yaxislimit = 250;

for i=1:k-1
    es(i) = hist(i).Es(4,2); % buyer idx, seller idx
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
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

leg1 = legend([{"$e^s_{24}$"}, {"$e^b_{24}$"}, {"$\bar{e}_{24}$"}, {"$z_{24}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
xticks([0 50 100 150])
yticks([0 50 100 150 200 250])
xlabel("Iterations")
ylabel('Energy (kWh)')
set(gca,'FontSize',15)

subplot(2,2,2)
pp1 = plot(1:(k-1),sellprice,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'LineWidth',1.5, 'color',[1 0.5 0]);
hold off
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

leg2 = legend([{"$\lambda^s_{24}$"}, {"$\lambda^b_{24}$"}]);
set(leg2,'Interpreter','latex','FontSize',16);
xlim([0 NumIterations])
ylim([0 8])
xticks([0 50 100 150])
yticks([0 2 4 6 8])
xlabel("Iterations")
ylabel("Price ("+char(0162)+"/kWh)")
set(gca,'FontSize',15)

for i=1:k-1
    es(i) = hist(i).Es(4,4);
    eb(i) = hist(i).Eb(4,4);
    eaverage(i) = hist(i).Eaverage(4,4);
    zglob(i) = hist(i).zglob(4,4);
    sellprice(i) = hist(i).Lambda_s(4,4);
    buyprice(i) = hist(i).Lambda_b(4,4);
end
subplot(2,2,3)
pp1 = plot(1:(k-1),es,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),eb,'-','LineWidth',1.5, 'color',[1 0.5 0]);
pp3 = plot(1:(k-1),eaverage,'-','LineWidth',1.5);
pp4 = plot(1:(k-1),zglob,'--','LineWidth',1.5,'color','k');
hold off
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

leg1 = legend([{"$e^s_{44}$"}, {"$e^b_{44}$"}, {"$\bar{e}_{44}$"}, {"$z_{44}$"}]);
set(leg1,'Interpreter','latex','FontSize',16,'numcolumns',2);
ylim([0 yaxislimit])
xlim([0 NumIterations])
yticks([0 50 100 150 200 250])
xticks([0 50 100 150])
ylabel('Energy (kWh)')
xlabel("Iterations")
set(gca,'FontSize',15)

subplot(2,2,4)
pp1 = plot(1:(k-1),sellprice,'LineWidth',1.5, 'color',[0 0.4 0.8]);
hold on
pp2 = plot(1:(k-1),buyprice,'LineWidth',1.5, 'color',[1 0.5 0]);
hold off
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)
leg2 = legend([{"$\lambda^s_{44}$"}, {"$\lambda^b_{44}$"}]);
set(leg2,'Interpreter','latex','FontSize',16);
xlim([0 NumIterations])
ylim([0 8])
xticks([0 50 100 150])
yticks([0 2 4 6 8])
ylabel("Price ("+char(0162)+"/kWh)")
xlabel("Iterations")
set(gca,'FontSize',15)

end