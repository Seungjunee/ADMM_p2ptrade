function [outputArg1,outputArg2] = main_result_graph(T, no, br, agents, sellers, buyers,...
                                    ADMM, ADMM_unreg, ADMM_pref, centralized_result, mpc, const, preference)
f1 = figure(1);
plot(1:br,abs(ADMM.Pline),'-','LineWidth',2,'color',[0 0.2 0.5]);
hold on
plot(1:br,abs(ADMM_unreg.Pline),'--','LineWidth',2,'color',[0 0.2 0.5]);
stairs(1:br,const.Linelimit,'--','LineWidth',2,'color',[1 0 0]);
hold off
xlabel('Line Index')
ylabel('Line power (MW)')
legend({'With Regularization','Without Regularization','Line Power Limits'})
ylim([0 5])
xlim([0 32])
xticks(0:4:32)
yticks([0 1 2 3 4 5])
set(gca,'FontSize',20)
set(f1,'Position',[0 0 1000 350])
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

f2 = figure(2);
plot(0:no-1,ADMM.VM,'-','LineWidth',2,'color',[0 0 0]);
hold on
plot(0:no-1,ADMM_unreg.VM,'--','LineWidth',2,'color',[0 0 0]);
plot(0:no-1,const.Vmin*ones(no,1),'--','LineWidth',2,'color',[1 0 0]);
plot(0:no-1,const.Vmax*ones(no,1),'--','LineWidth',2,'color',[1 0 0]);
hold off
xlabel('Node Index')
ylabel('Voltage (p.u.)')
legend({'With Regularization','Without Regularization','Voltage Limits'})
xlim([0 32])
ylim([0.92 1.08])
xticks(0:4:32)
yticks([0.92 0.96 1.00 1.04 1.08])
% yticks([0.95 1.00 1.05])
set(gca,'FontSize',20)
set(f2,'Position',[0 0 1000 350])
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

figure(3)
bar(1:agents.no,[[centralized_result.selltrade; centralized_result.buytrade], [ADMM.sell; ADMM.buy]]);
xlabel('agent index','FontSize',12)
ylabel('Net trade power [kW]','FontSize',12)
legend({'centralized','ADMM w/ reg.'})
set(gca,'FontSize',16)

figure(4) % except retailer
buyer_idx = 4;
% target_idx = find(ADMM_pref.energy(1:end,agent_idx)>1e-3)';
target_idx = buyers(buyer_idx).partner;
aH = axes;
aH.GridAlpha=1;
aH.GridColor=[0.8 0.8 0.8];
system_fee = [zeros(length(target_idx),2),min(ADMM_pref.sellprice(buyer_idx,target_idx)',ADMM_pref.buyprice(buyer_idx,target_idx)')];
maximum = [ADMM_pref.sellprice(buyer_idx,target_idx)', ADMM_pref.buyprice(buyer_idx,target_idx)',...
           max(ADMM_pref.sellprice(buyer_idx,target_idx)',ADMM_pref.buyprice(buyer_idx,target_idx)')];
bH = bar(maximum);
bH(1).FaceColor = [0 0.4 0.8];
bH(2).FaceColor = [1 0.5 0];
bH(3).FaceColor = [1 0.23 0.23];
hold on;
bH1 = bar(system_fee);
bH1(3).BarWidth = 1.2;
for ii = bH1
    ii.FaceColor = [1 1 1];
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
for ii = bH
    ii.LineStyle = 'None';
    ii.EdgeColor = [1 1 1];
end
p4 = scatter((1:length(target_idx))+0.225,(ADMM_pref.sellprice(buyer_idx,target_idx)'+ADMM_pref.buyprice(buyer_idx,target_idx)')/2,100,'k*');
hold off
xlabel('seller index','FontSize',12)
ylabel("Buyer 4's price ("+char(0162)+"/kWh)",'FontSize',12)
leg1 = legend([bH,p4],{'$\lambda^s_{ij}$','$\lambda^b_{ij}$','2$\tilde{\lambda}_{ij}$','$\bar{\lambda}_{ij}$'});
set(leg1,'Interpreter','latex','FontSize',24)
set(gca,'FontSize',24)
set(gca,'xticklabel',target_idx)
ylim([0 10])
yticks([0 2 4 6 8 10])
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

figure(45) % except retailer
buyer_idx = 4;
% target_idx = find(ADMM_pref.energy(1:end,agent_idx)>1e-3)';
target_idx = buyers(buyer_idx).partner;
aH = axes;
aH.GridAlpha=1;
aH.GridColor=[0.8 0.8 0.8];
system_fee = [zeros(length(target_idx),2),ADMM_pref.buyprice(buyer_idx,target_idx)'-preference(buyer_idx,target_idx)'];
maximum = [ADMM_pref.sellprice(buyer_idx,target_idx)', ADMM_pref.buyprice(buyer_idx,target_idx)',...
           ADMM_pref.buyprice(buyer_idx,target_idx)'];
bH = bar(maximum);
bH(1).FaceColor = [0 0.4 0.8];
bH(2).FaceColor = [1 0.5 0];
bH(3).FaceColor = [1 0.23 0.23];
hold on;
bH1 = bar(system_fee);
bH1(3).BarWidth = 1.2;
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
ylim([0 10])
leg1 = legend([bH],{'$\lambda^s_{ij}$','$\lambda^b_{ij}$','$u_{ij}$'});
set(leg1,'Interpreter','latex','FontSize',24)
set(gca,'FontSize',24)
set(gca,'xticklabel',target_idx)
yticks([0 2 4 6 8 10])
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)



figure(44) % except retailer
buyer_idx = 4;
% target_idx = find(ADMM_pref.energy(1:end,agent_idx)>1e-3)';
target_idx = buyers(buyer_idx).partner;
aH = axes;
homogeneous_price = [zeros(length(target_idx),1),ADMM_pref.buyprice(buyer_idx,target_idx)'-preference(buyer_idx,target_idx)'];
maximum = [ADMM_pref.sellprice(buyer_idx,target_idx)', ADMM_pref.buyprice(buyer_idx,target_idx)'];
plot([0:5],ones(1,6)*mean(ADMM_pref.buyprice(buyer_idx,target_idx)-preference(buyer_idx,target_idx)),'--k')
hold on;
bH = bar(target_idx,maximum);
bH(1).BarWidth = 0.6;
bH(2).BarWidth = 0.6;
bH(1).FaceColor = [0 0.4 0.8];
bH(2).FaceColor = [1 0.75 0.05]; % 0.2 0.8 0.8
bH2 = bar(target_idx,homogeneous_price);
bH2(1).BarWidth = 0.6;
bH2(1).BarWidth = 0.6;
for ii = bH2
    ii.FaceColor = [0.95 0.25 0.05];
    ii.EdgeColor = [0 0 0];
end
for ii = bH
    ii.EdgeColor = [0 0 0];
end
hold off
xlabel('seller index','FontSize',12)
ylabel("Buyer 4's price ("+char(0162)+"/kWh)",'FontSize',12)
leg1 = legend([bH, bH2(2)],{'$\lambda^s_{ij}$', "$u_{ij}$",'$\lambda^b_{ij}-u_{ij}$'});
set(leg1,'Interpreter','latex','FontSize',24)
set(gca,'FontSize',24)
xticks(target_idx)
ylim([0 10])
yticks([0 2 4 6 8 10])
xlim([1.5 4.5])
grid on
set(gca,'Gridcolor',[0.9 0.9 0.9],'GridAlpha',1)

figure(5) % net trade power seller & buyer
subplot(1,2,1)
bar(1:length(sellers),[centralized_result.selltrade,ADMM.sell]);
xlabel('seller index')
ylabel('Trade energy[kW]')
ylim([0 200])
legend({'ADMM','centralized'})
set(gca,'FontSize',16)

subplot(1,2,2)
bar(1:length(buyers),[centralized_result.buytrade,ADMM.buy]);
xlabel('buyer index')
legend({'ADMM','centralized'})
ylim([0 200])
set(gca,'FontSize',16)

f6 = figure(6);

ADMM.energy(ADMM.energy<1e-2) = 0; % invisible line for less than 10W
grid_price = ADMM.buyprice-ADMM.sellprice;
p=plot(T,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'NodeColor','k','EdgeColor','k','LineWidth',1,'EdgeAlpha',1);
hold on
diG = digraph; 
G1 = addnode(diG,T.Nodes); % normal
G2 = addnode(diG,T.Nodes); % positive
G3 = addnode(diG,T.Nodes); % negative
energy_graph = [];
physics_energy = blkdiag(0,agents.Ab*ADMM.energy*agents.As');
physics_price = blkdiag(0,agents.Ab*grid_price*agents.As');
for j=1:no
    for i=1:no
        if physics_energy(j,i)>=10
            if physics_price(j,i) > 1e-3
                G2 = addedge(G2,i, j,physics_energy(j,i));
            elseif physics_price(j,i) < -1e-3
                G3 = addedge(G3,i, j,physics_energy(j,i));
            else
                G1 = addedge(G1,i, j,physics_energy(j,i));
            end
        end
    end
end
if ~isempty(G1.Edges)
G1LWidths =G1.Edges.Weight/12;
p1 = plot(G1,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[0.5, 0.5, 0.5],'LineWidth',G1LWidths,'EdgeAlpha',0.5);
p1.ArrowSize = 10;
end
if ~isempty(G2.Edges)
G2LWidths =G2.Edges.Weight/12;
p2 = plot(G2,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[0.3, 0.75, 1],'LineWidth',G2LWidths,'EdgeAlpha',0.5);
p2.ArrowSize = 10;
end
if ~isempty(G3.Edges)
G3LWidths =G3.Edges.Weight/12;
p3 = plot(G3,'XData',mpc.XY_data(:,2),'YData',mpc.XY_data(:,3),'EdgeColor',[1, 0.75, 0],'LineWidth',G3LWidths,'EdgeAlpha',0.5);
% ,'EdgeAlpha',1
p3.ArrowSize = 10;
end
p3.NodeColor = 'k';
% highlight(p3,[sellers.bus],'NodeColor',[70 114 196]/255)
% highlight(p3,[buyers.bus],'NodeColor',[237 125 49]/255)

p.NodeLabel = {};
p1.NodeLabel = {};
p2.NodeLabel = {};
p3.NodeLabel = {};
hold off
set(gca, 'Ydir', 'reverse')
set(gca,'xtick',[])
set(gca,'ytick',[])
set(f6,'Position',[0 0 1000 250])

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

