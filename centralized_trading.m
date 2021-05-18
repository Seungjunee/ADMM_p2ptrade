function central = centralized_trading(agents, sellers, buyers, mpc, const)
%% Network data
no = length(mpc.bus);
br = length(mpc.branch);

%% Load data read

Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Ybus
Zbus = inv(Ybus(2:no,2:no))';

%% Optimal Power Flow problem
%--------------------------------------------------------------------------
incidence = zeros(length(sellers),length(buyers));
for i = 1:length(sellers)
    incidence(i,sellers(i).partner)=1;
end

[f, qf, A, b, Aeq, beq, ub, lb] = deal([]); % initialization

pos.sell_trade = [0,size(sellers,2)];
pos.buy_trade = [pos.sell_trade(end),pos.sell_trade(end)+size(buyers,2)];
pos.bilateral_trade = [pos.buy_trade(end),pos.buy_trade(end)+length(sellers)*length(buyers)];
dl = pos.bilateral_trade(end);

% quad cost of agents
qf = zeros(1,dl);
for i=1:size(sellers,2)
    qf(pos.sell_trade(1)+i) = sellers(i).coeff.A;
end
for j=1:size(buyers,2)
    qf(pos.buy_trade(1)+j) = buyers(j).coeff.A;
end

% linear cost of agents
f = zeros(1,dl);
for i=1:size(sellers,2)
    f(pos.sell_trade(1)+i) = sellers(i).coeff.B;
end
for j=1:size(buyers,2)
    f(pos.buy_trade(1)+j) = -buyers(j).coeff.B;
end
% for pool based balance constraint

for i=1:size(sellers,2)
    ev = zeros(1,dl);
    temp = incidence;
    temp(setdiff(1:length(sellers),i),:) = 0;
    ev(pos.sell_trade(1)+i) = 1;
    ev(pos.bilateral_trade(1)+find(temp')) = -1;
    Aeq = [Aeq; ev];
    beq = [beq; 0];
end
for j=1:size(buyers,2)
    ev = zeros(1,dl);
    temp = incidence';
    temp(setdiff(1:length(buyers),j),:) = 0;
    ev(pos.buy_trade(1)+j) = 1;
    ev(pos.bilateral_trade(1)+find(temp)) = -1;
    Aeq = [Aeq; ev];
    beq = [beq; 0];
end


% sell_trade sell_util sell_net buy_trade buy_util buy_net

lb = [];
ub = [];

for i=1:size(sellers,2)
    lb = [lb, sellers(i).gen];
    ub = [ub, sellers(i).genmax];
end

for j=1:size(buyers,2)
    lb = [lb, buyers(j).load];
    ub = [ub, buyers(j).loadmax];
end
lb = [lb, zeros(1, size(sellers,2)*size(buyers,2))];
ub = [ub, 1000*reshape(incidence',1,[])];
%% Network constraints
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));
Pline_fix = pfresult.branch(:,14);
Qline_fix = pfresult.branch(:,15);
V_fix = pfresult.bus(:,8);
ISF=makePTDF(mpc);

%% Line loading constraint
if const.activate_Linelimit == true
    Lineconst = zeros(br,dl);
    for i=1:length(sellers)
        Lineconst(1:br,pos.sell_trade(1)+i) = +ISF(:,sellers(i).bus);
    end
    for j=1:length(buyers)
        Lineconst(1:br,pos.buy_trade(1)+j) = -ISF(:,buyers(j).bus);
    end
    Lineconst = Lineconst/1e3;
    A = [A; Lineconst;-Lineconst];
    b = [b; (const.Linelimit-Pline_fix);(const.Linelimit+Pline_fix)];
end
%% Voltage constraint
if const.activate_Voltagelimit == true
    Voltageconst = zeros(no-1,dl);
    for i=1:length(sellers)
        if sellers(i).bus == 1 % for slack bus
            Voltageconst(1:no-1,pos.sell_trade(1)+i) =  zeros(no-1,1);
            continue
        end
        Voltageconst(1:no-1,pos.sell_trade(1)+i) = -real(Zbus(:,sellers(i).bus-1))/mpc.baseMVA;
    end
    for j=1:length(buyers)
        if buyers(j).bus == 1 % for slack bus
            Voltageconst(1:no-1,pos.buy_trade(1)+j) =  zeros(no-1,1);
            continue
        end
        Voltageconst(1:no-1,pos.buy_trade(1)+j) = +real(Zbus(:,buyers(j).bus-1))/mpc.baseMVA;
    end
    Voltageconst = Voltageconst/1e3;
    A = [A; Voltageconst; -Voltageconst];
    b = [b; V_fix(2:end)-const.Vmin*ones(no-1,1);const.Vmax*ones(no-1,1)-V_fix(2:end)];
end
%% OPF process...
[x,fval,exitflag,output,lambda] = quadprog(diag(qf)*2,f,A,b,Aeq,beq,lb,ub); % quadratic term 1/2 in process
central.selltrade = x(pos.sell_trade(1)+1:pos.sell_trade(end));
central.buytrade = x(pos.buy_trade(1)+1:pos.buy_trade(end));
central.sellprice = lambda.eqlin;
central.buyprice = lambda.eqlin;

Trade_map = agents.As*central.selltrade-agents.Ab*central.buytrade;
Trade_map = Trade_map/1000;
