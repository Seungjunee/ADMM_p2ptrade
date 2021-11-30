function central = centralized_trading(agents, sellers, buyers, mpc, const)
%% Network data
no = length(mpc.bus);
br = length(mpc.branch);

%% Load data read

Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Ybus
Zbus = inv(Ybus(2:no,2:no))';

pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));
VM = pfresult.bus(2:end,8); % except slack bus
VA = pfresult.bus(2:end,9)*pi/180;
alpha = zeros(no-1,no-1);
beta = zeros(no-1,no-1);
for i = 1:no-1
    for j = 1:no-1
        alpha(i,j) = real(Zbus(i,j))/(VM(i)*VM(j))*cos(VA(i)-VA(j));
        beta(i,j) = real(Zbus(i,j))/(VM(i)*VM(j))*sin(VA(i)-VA(j));
    end
end
LSF = 2*((alpha*(-pfresult.bus(2:end,3))-beta*(-pfresult.bus(2:end,4))))/mpc.baseMVA;
MLSFji = LSF'*agents.As - (LSF'*agents.Ab)';

%% Optimal Power Flow problem
%--------------------------------------------------------------------------
incidence = zeros(length(sellers),length(buyers));
for i = 1:length(sellers)
    incidence(i,sellers(i).partner)=1;
end
cij_vec = MLSFji(:);
cij_vec(cij_vec>0) = 7 * cij_vec(cij_vec>0);
cij_vec(cij_vec<0) = 3 * cij_vec(cij_vec<0);


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

if const.activate_Loss == true
for i=1:size(sellers,2)*size(buyers,2)
    f(pos.bilateral_trade(1)+i) = f(pos.bilateral_trade(1)+i)+cij_vec(i);
end
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
%% OPF process with Gurobi
model.Q = sparse(diag(qf));
model.A = sparse([A;Aeq]);
model.rhs = [b;beq];
model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
params.outputflag = 0;
model.obj = f;
model.lb = lb;
model.ub = ub;

results = gurobi(model, params);
% [x,fval,exitflag,output,lambda] = quadprog(diag(qf)*2,f,A,b,Aeq,beq,lb,ub); % quadratic term 1/2 in process
x=results.x;
central.selltrade = x(pos.sell_trade(1)+1:pos.sell_trade(end));
central.buytrade = x(pos.buy_trade(1)+1:pos.buy_trade(end));
central.sellprice = results.pi(size(A,1)+(1:length(sellers)));
central.buyprice = results.pi(size(A,1)+length(sellers)+1:end);

Trade_map = agents.As*central.selltrade-agents.Ab*central.buytrade;
Trade_map = Trade_map/1000;
