function [outputArg1,outputArg2] = main_result_graph(T, no, br, agents, sellers, buyers,...
                                    ADMM, ADMM_unreg, ADMM_pref, centralized_result, mpc, const, preference)
f1 = figure(1);
plot(1:br,abs(ADMM.Pline),'-','LineWidth',2,'color',[0 0.2 0.5]);
hold on
plot(1:br,abs(ADMM_unreg.Pline),'--','LineWidth',2,'color',[0 0.2 0.5]);
stairs(1:br,const.Linelimit,'--','LineWidth',2,'color',[1 0 0]);
hold off
xlabel('branch index')
ylabel('Line power(MW)')
legend({'w/ regulation','w/o regulation','line power limit'})
ylim([0 2.5])
xlim([1 br])
set(gca,'FontSize',24)
set(f1,'Position',[0 0 1000 500])
grid on

f2 = figure(2);
plot(1:no,ADMM.VM,'-','LineWidth',2,'color',[0 0 0]);
hold on
plot(1:no,ADMM_unreg.VM,'--','LineWidth',2,'color',[0 0 0]);
plot(1:no,const.Vmin*ones(no,1),'--','LineWidth',2,'color',[1 0 0]);
hold off
xlabel('bus index')
ylabel('Voltage magnitude (p.u.)')
legend({'w/ regulation','w/o regulation','voltage limit'})
xlim([1 no])
ylim([0.94 1.01])
yticks([0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01])
set(gca,'FontSize',24)
set(f2,'Position',[0 0 1000 500])
grid on

figure(3)
bar(1:agents.no,[[centralized_result.selltrade; centralized_result.buytrade], [ADMM.sell; ADMM.buy]]);
xlabel('agent index','FontSize',12)
ylabel('Net trade power [kW]','FontSize',12)
legend({'centralized','ADMM w/ reg.'})
ylim([0 200])
set(gca,'FontSize',16)

figure(4) % except retailer
agent_idx = 4;
% target_idx = find(ADMM_pref.energy(1:end,agent_idx)>1e-3)';
target_idx = buyers(agent_idx).partner;
aH = axes;
system_fee = [zeros(length(target_idx),2),min(ADMM_pref.sellprice(target_idx,agent_idx),ADMM_pref.buyprice(target_idx,agent_idx))];
homogeneous_price = [zeros(length(target_idx),1),...
    ADMM_pref.buyprice(target_idx,agent_idx)-preference(target_idx,agent_idx),zeros(length(target_idx),1)];
maximum = [ADMM_pref.sellprice(target_idx,agent_idx), ADMM_pref.buyprice(target_idx,agent_idx),...
           max(ADMM_pref.sellprice(target_idx,agent_idx),ADMM_pref.buyprice(target_idx,agent_idx))];
bH = bar(maximum);
bH(1).FaceColor = [0 0.4 0.8];
bH(2).FaceColor = [1 0.5 0];
bH(3).FaceColor = [1 0.23 0.23];
hold on;
bH1 = bar(system_fee);
bH2 = bar(homogeneous_price);
% bH1(3).Visible = 'off';
for ii = bH2
    ii.FaceColor = [0.6 0.45 0.2];
    ii.LineStyle = 'None';
    ii.EdgeColor = [0.6 0.45 0.2];
end
for ii = bH1
    ii.FaceColor = [1 1 1];
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
for ii = bH
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
hold off
xlabel('seller index','FontSize',12)
ylabel("Buyer 4's price ("+char(0162)+"/kWh)",'FontSize',12)
leg1 = legend([bH, bH2(2)],{'$\lambda^s_{ij}$', "$c_{ij}$", '2$\bar{\lambda}_{ij}$', '$\mu^b_{j}$'});
set(leg1,'Interpreter','latex','FontSize',16)
set(gca,'FontSize',20)
set(gca,'xticklabel',target_idx-1)
ylim([0 10])
grid on

figure(44) % except retailer
agent_idx = 3+1;
target_idx = buyers(agent_idx).partner;
aH = axes;
system_fee = [zeros(length(target_idx),2),min(ADMM_pref.sellprice(target_idx,agent_idx),ADMM_pref.buyprice(target_idx,agent_idx))];
homogeneous_price = [zeros(length(target_idx),1),...
    ADMM_pref.buyprice(target_idx,agent_idx)-preference(target_idx,agent_idx),zeros(length(target_idx),1)];
maximum = [ADMM_pref.sellprice(target_idx,agent_idx), ADMM_pref.buyprice(target_idx,agent_idx),...
           max(ADMM_pref.sellprice(target_idx,agent_idx),ADMM_pref.buyprice(target_idx,agent_idx))];
bH = bar(maximum);
bH(1).FaceColor = [0 0.4 0.8];
bH(2).FaceColor = [1 0.5 0];
bH(3).FaceColor = [1 0.23 0.23];
hold on;
bH1 = bar(system_fee);
bH2 = bar(homogeneous_price);
% bH1(3).Visible = 'off';
for ii = bH2
    ii.FaceColor = [0.6 0.45 0.2];
    ii.LineStyle = 'None';
    ii.EdgeColor = [0.6 0.45 0.2];
end
for ii = bH1
    ii.FaceColor = [1 1 1];
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
for ii = bH
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
hold off
xlabel('seller index','FontSize',12)
ylabel("Buyer 3's price ("+char(0162)+"/kWh)",'FontSize',12)
leg1 = legend([bH, bH2(2)],{'$\lambda^s_{ij}$', "$c_{ij}$", '2$\bar{\lambda}_{ij}$', '$\mu^b_{j}$'});
set(leg1,'Interpreter','latex','FontSize',16)
set(gca,'FontSize',20)
set(gca,'xticklabel',target_idx-1)
ylim([0 10])
grid on


figure(5) % net trade power seller & buyer
subplot(1,2,1)
bar(1:length(sellers)-1,[centralized_result.selltrade(2:end),ADMM.sell(2:end)]);
xlabel('seller index')
ylabel('Trade energy[kW]')
ylim([0 200])
legend({'ADMM','centralized'})
set(gca,'FontSize',16)

subplot(1,2,2)
bar(1:length(buyers)-1,[centralized_result.buytrade(2:end),ADMM.buy(2:end)]);
xlabel('buyer index')
legend({'ADMM','centralized'})
ylim([0 200])
set(gca,'FontSize',16)

f6 = figure(6);

ADMM.energy(ADMM.energy<1e-2) = 0; % invisible line for less than 10W
grid_price = ADMM.buyprice-ADMM.sellprice;
p=plot(T,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'NodeColor','k','EdgeColor','k','LineWidth',1);
hold on
diG = digraph; 
G1 = addnode(diG,T.Nodes); % normal
G2 = addnode(diG,T.Nodes); % positive
G3 = addnode(diG,T.Nodes); % negative
energy_graph = [];
physics_energy = blkdiag(0,agents.As*ADMM.energy*agents.Ab');
physics_price = blkdiag(0,agents.As*grid_price*agents.Ab');
for i=1:no
    for j=1:no
        if physics_energy(i,j)~=0
            if physics_price(i,j) > 1e-3
                G2 = addedge(G2,i, j,physics_energy(i,j));
            elseif physics_price(i,j) < -1e-3
                G3 = addedge(G3,i, j,physics_energy(i,j));
            else
                G1 = addedge(G1,i, j,physics_energy(i,j));
            end
        end
    end
end
if ~isempty(G1.Edges)
G1LWidths =G1.Edges.Weight/8;
p1 = plot(G1,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[0.5, 0.5, 0.5],'LineWidth',G1LWidths);
p1.ArrowSize = 10;
end
if ~isempty(G2.Edges)
G2LWidths =G2.Edges.Weight/8;
p2 = plot(G2,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[0.07, 0.62, 1],'LineWidth',G2LWidths);
p2.ArrowSize = 10;
end
if ~isempty(G3.Edges)
G3LWidths =G3.Edges.Weight/8;
p3 = plot(G3,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[1, 0, 0],'LineWidth',G3LWidths);
p3.ArrowSize = 10;
end
p3.NodeColor = 'k';
p.NodeLabel = {};
p1.NodeLabel = {};
p2.NodeLabel = {};
p3.NodeLabel = {};
hold off
set(gca, 'Ydir', 'reverse')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(f6,'Position',[0 0 1000 300])

% f7 = figure(7);
% 
% mpc.bus(:,3) = mpc.bus(:,3)-sum(physics_energy,2)/1e3;
% mpc.bus(:,3) = mpc.bus(:,3)+sum(physics_energy,1)'/1e3;
% pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));
% 
% plot(1:br,abs(ADMM.Pline),'-','LineWidth',2,'color',[0 0.2 0.5]);
% hold on
% plot(1:br,abs(pfresult.branch(:,14)),'--','LineWidth',2,'color',[0 0.2 0.5]);
% stairs(1:br,const.Linelimit,'--','LineWidth',2,'color',[1 0 0]);
% hold off
% xlabel('branch index')
% ylabel('Line power(MW)')
% legend({'w/ regulation','w/o regulation','line power limit'})
% ylim([0 2.5])
% xlim([1 br])
% set(gca,'FontSize',24)
% set(f1,'Position',[0 0 1000 500])
% grid on
% 
% f8 = figure(8);
% plot(1:no,ADMM.VM,'-','LineWidth',2,'color',[0 0 0]);
% hold on
% plot(1:no,pfresult.bus(:,8),'--','LineWidth',2,'color',[0 0 0]);
% plot(1:no,const.Vmin*ones(no,1),'--','LineWidth',2,'color',[1 0 0]);
% hold off
% xlabel('bus index')
% ylabel('Voltage magnitude (p.u.)')
% legend({'w/ regulation','w/o regulation','voltage limit'})
% xlim([1 no])
% ylim([0.94 1.01])
% yticks([0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01])
% set(gca,'FontSize',24)
% set(f2,'Position',[0 0 1000 500])
% grid on
end

