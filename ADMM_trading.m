function [Eaverage, Lambda_s, Lambda_b] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference)
%% P2P Energy Trading algorithm with consensus ADMM
%--------------------------------------------------------------------------
% This algorithm is implemented with 'Gurobi' optimization solver and
% 'matpower' package
% Test feeder is modified IEEE 33bw
%--------------------------------------------------------------------------
%% Network data

no = length(mpc.bus);  % the number of bus
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));

%% Load data read

Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Ybus
Zbus = inv(Ybus(2:no,2:no))';

%% Energy trading process...
% initiatlization
Es = zeros(length(sellers),length(buyers));
Eb = zeros(length(sellers),length(buyers));
zglob = zeros(length(sellers),length(buyers));
Eaverage = zeros(length(sellers),length(buyers));
% price is initialized by feed-in-tariffs
Lambda_b = buyers(1).coeff.B*ones(length(sellers),length(buyers));
Lambda_s = buyers(1).coeff.B*ones(length(sellers),length(buyers));

%% market PTDFs
MarketPTDF = [];

if const.activate_Linelimit == true
    ISF = makePTDF(mpc,1); % matpower package support injection shift factor
    ISF(:,1) = []; % reduced first column corresponding slack bus
    sell_ISF = ISF*agents.As;
    buy_ISF = ISF*agents.Ab;
    for j=1:length(buyers)
        for i=1:length(sellers)
%             MarketPTDF = [MarketPTDF, ISF(:,sellers(i).bus)-ISF(:,buyers(j).bus)];
            MarketPTDF = [MarketPTDF, sell_ISF(:,i)-buy_ISF(:,j)];
        end
    end
    MarketPTDF(abs(MarketPTDF)<1e-3)=0;
end
MarketPTDF = MarketPTDF/1e3; % convert kW to MW
%% market MVSCs
MarketMVSC = [];
if const.activate_Voltagelimit == true
    sell_VSF = real(Zbus)*agents.As;
    buy_VSF = real(Zbus)*agents.Ab;
    for j=1:length(buyers)
        for i=1:length(sellers)
            MarketMVSC = [MarketMVSC, (sell_VSF(:,i)-buy_VSF(:,j))/mpc.baseMVA];
        end
    end
end
MarketMVSC = MarketMVSC/1e3; % convert kW to MW
%% market MLSFs
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

MLSFij = (LSF'*agents.As)' - LSF'*agents.Ab;
for j=1:length(buyers)
    for i=1:length(sellers)
        if sellers(i).bus == 1 || buyers(j).bus == 1 % for slack bus / neither seller/buyer
            MLSFij(i,j) = 0; % We impose loss transaction fee on retailer..
            continue
        end
%         MLSFij(i,j) = (LSF(sellers(i).bus)-LSF(buyers(j).bus));
        if MLSFij(i,j) > 0 % 
            cij(i,j) = sellers(1).coeff.B*MLSFij(i,j);
        else
            cij(i,j) = -buyers(1).coeff.B*MLSFij(i,j);
        end
    end
end
cij_vec = reshape(cij,[],1);
if const.activate_Loss == false % if not activate loss, set to zeros
    cij_vec = zeros(size(cij_vec));
end
%% For fixed load(nonparticipate p2p energy trade)
%--------------------------------------------------------------------------
Pline_fix = pfresult.branch(:,14);
Qline_fix = pfresult.branch(:,15);
V_fix = pfresult.bus(:,8);
%--------------------------------------------------------------------------
% DSO's optimization model (safe trading set) -----------------------------
model.Q = sparse(eye(length(sellers)*length(buyers)));
model.A = sparse([MarketPTDF; -MarketPTDF; -MarketMVSC]);
model.rhs = [];
if const.activate_Linelimit == true
    model.rhs = [(const.Linelimit-Pline_fix); (const.Linelimit+Pline_fix)];
end
if const.activate_Voltagelimit == true
    model.rhs = [model.rhs; V_fix(2:end)-const.Vmin*ones(no-1,1)];
end
model.sense = repmat('<',size(model.A,1),1);
params.outputflag = 0;
%--------------------------------------------------------------------------
for round=1:1
    k=1;    % iteration index
    epsilon_pri = 1;    % primal residual
    epsilon_dual = 1;   % dual residual
    zglob_prev = 0;     % initialize previous global variable values
    
    while (epsilon_pri > 1e-5) || (epsilon_dual > 1e-3)
        % local variable (amount) update for sellers
        for i=1:length(sellers)
            lambdas = Lambda_s(i,sellers(i).partner);
            eglobs = zglob(i,sellers(i).partner);
            Es(i,sellers(i).partner) = sellers(i).local_opt(eglobs,lambdas);
        end
        % local variable (amount) update for buyers
        for j=1:length(buyers)
            lambdab = Lambda_b(buyers(j).partner,j)';
            eglobb = zglob(buyers(j).partner,j)';
            preference_j = preference(buyers(j).partner,j)';
            Eb(buyers(j).partner,j) = buyers(j).local_opt(eglobb,lambdab,preference_j)';
        end
        % average energy amount
        Eaverage = (Es+Eb)./2;
        Eaverage_vec = reshape(Eaverage,[],1);
        % difference of energy price
        Lambda_bar = (Lambda_b-Lambda_s)./2;
        Lambda_barvec = reshape(Lambda_bar,[],1);
        
        % global energy variable (amount) update
        if const.activate_Linelimit == true || const.activate_Voltagelimit == true
            model.obj = -2*Eaverage_vec+(-2*Lambda_barvec+cij_vec)./rho;
            results = gurobi(model, params);
            z = results.x;
            zglob = reshape(z,length(sellers),length(buyers));
        elseif const.activate_Loss == true
            z = Eaverage_vec+(Lambda_barvec-cij_vec/2)./rho;
            zglob = reshape(z,length(sellers),length(buyers));
        else
            zglob = Eaverage;
        end
        
        % dual variable (price) update
        Lambda_s = Lambda_s + rho*(zglob-Es);
        Lambda_b = Lambda_b + rho*(Eb-zglob);
        
        % residuals update(stopping criteria)
        epsilon_pri = norm(Eaverage-zglob,2);
        epsilon_dual = rho*norm(zglob-zglob_prev,2);
        
        zglob_prev = zglob;
        hist(k).Es = Es;
        hist(k).Eb = Eb;
        hist(k).zglob = zglob;
        hist(k).Eaverage = Eaverage;
        hist(k).Lambda_s = Lambda_s;
        hist(k).Lambda_b = Lambda_b;
        hist(k).primal = epsilon_pri;
        hist(k).dual = epsilon_dual;
        k=k+1;
    end
%     ADMM_graph_result(hist, k, MLSFij)
    % derive utility of overall market players
    utility = 0;
    for j=1:length(buyers)
        utility = utility + (buyers(j).coeff.B*sum(Eb(:,j))-buyers(j).coeff.A*sum(Eb(:,j))^2);
    end
    for i=1:length(sellers)
        utility = utility - (sellers(i).coeff.B*sum(Es(i,:))+sellers(i).coeff.A*sum(Es(i,:))^2);
    end
    pi_loss = sum(MLSFij.*(Es),'all');
    pi_loss_fees = sum(cij.*(Es),'all');
    system_fees = sum((Lambda_b-Lambda_s).*(Es),'all');
end
end

