function [Eaverage, Lambda_s, Lambda_b] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference)
%% P2P Energy Trading algorithm with consensus ADMM
%--------------------------------------------------------------------------
% This algorithm is implemented with 'Gurobi' optimization solver and
% 'matpower' package
% Test feeder is modified IEEE 33bw
%--------------------------------------------------------------------------
%% Network data

no = length(mpc.bus);  % the number of bus
br = length(mpc.branch);
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));

incidence = zeros(length(buyers),length(sellers));
for i = 1:length(sellers)
    incidence(sellers(i).partner,i)=1;
end
tradelen = sum(incidence,'all');
incidence = logical(incidence);
%% Load data read

Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Ybus
Zbus = inv(Ybus(2:no,2:no))';

%% Energy trading process...
% initiatlization
Es = zeros(length(buyers),length(sellers));
Eb = zeros(length(buyers),length(sellers));
zglob = zeros(length(buyers),length(sellers));
Eaverage = zeros(length(buyers),length(sellers));
Lambda_b = zeros(length(buyers),length(sellers));
Lambda_s = zeros(length(buyers),length(sellers));

fd_low = zeros(length(br),1);
fd_upp = zeros(length(br),1);
vd_low = zeros(length(no-1),1);
vd_upp = zeros(length(no-1),1);
%% market PTDFs
MarketPTDF = [];
ISF = makePTDF(mpc,1); % matpower package support injection shift factor
ISF(:,1) = []; % reduced first column corresponding slack bus

if const.activate_Linelimit == true
    sell_ISF = ISF*agents.As;
    buy_ISF = ISF*agents.Ab;
    for i=1:length(sellers)
        for j=1:length(sellers(i).partner)
            MarketPTDF = [MarketPTDF, sell_ISF(:,i)-buy_ISF(:,sellers(i).partner(j))];
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
    for i=1:length(sellers)
        for j=1:length(sellers(i).partner)
            MarketMVSC = [MarketMVSC, (sell_VSF(:,i)-buy_VSF(:,sellers(i).partner(j)))/mpc.baseMVA];
        end
    end
    MarketPTDF(abs(MarketPTDF)<1e-3)=0;
end
MarketMVSC = MarketMVSC/1e3; % convert kW to MW
%% market MLSFs

LSF = zeros(no-1,1);
for i=1:no-1
    for l=1:32
        LSF(i) = LSF(i) + 2*pfresult.branch(l,3)*pfresult.branch(l,14)*sum(ISF(l,i))/100;
    end
end
MLSFji = LSF'*agents.As - (LSF'*agents.Ab)';
cij_vec = MLSFji(incidence);
cij_vec(cij_vec>0) = 7 * cij_vec(cij_vec>0); % retail price
cij_vec(cij_vec<0) = 3 * cij_vec(cij_vec<0); % fit price

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
model.Q = sparse(eye(tradelen));
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
    
    while (epsilon_pri > 1e-4) || (epsilon_dual > 1e-4)
        % local variable (amount) update for sellers
        for i=1:length(sellers)
            lambdas = Lambda_s(sellers(i).partner,i)';
            eglobs = zglob(sellers(i).partner,i)';
            Es(sellers(i).partner,i) = sellers(i).local_opt(eglobs,lambdas);
        end
        % local variable (amount) update for buyers
        for j=1:length(buyers)
            lambdab = Lambda_b(j,buyers(j).partner);
            eglobb = zglob(j,buyers(j).partner);
            preference_j = preference(j,buyers(j).partner);
            Eb(j,buyers(j).partner) = buyers(j).local_opt(eglobb,lambdab,preference_j);
        end
        % average energy amount
        Eaverage = (Es+Eb)./2;
        Eaverage_vec = Eaverage(incidence);
        % difference of energy price
        Lambda_bar = (Lambda_b-Lambda_s)./2;
        Lambda_barvec = Lambda_bar(incidence);
        
        % global energy variable (amount) update
        if const.activate_Linelimit == true || const.activate_Voltagelimit == true
            ehat = Eaverage_vec+(Lambda_barvec-cij_vec/2)./rho;
            model.obj = -2*Eaverage_vec+(-2*Lambda_barvec+cij_vec)./rho;
            results = gurobi(model, params);
            z = results.x;
            zglob(incidence) = z;
        elseif const.activate_Loss == true
            ehat = Eaverage_vec+(Lambda_barvec-cij_vec/2)./rho;
            zglob(incidence) = ehat;
        else
            zglob = Eaverage;
        end
        
        % dual variable (price) update
        Lambda_s = Lambda_s + rho*(zglob-Es);
        Lambda_b = Lambda_b + rho*(Eb-zglob);
        
        % residuals update(stopping criteria)
        epsilon_pri = norm([Es-zglob;Eb-zglob],2);
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
    ADMM_graph_result(hist, k)
    % derive utility of overall market players
    utility = 0;
    for j=1:length(buyers)
        utility = utility + (buyers(j).coeff.B*sum(Eb(j,:))-buyers(j).coeff.A*sum(Eb(j,:))^2);
    end
    for i=1:length(sellers)
        utility = utility - (sellers(i).coeff.B*sum(Es(:,i))+sellers(i).coeff.A*sum(Es(:,i))^2);
    end
    Ploss = sum(MLSFji.*(Eaverage),'all'); % kW
    pi_loss_fees = sum(cij_vec'*Eaverage_vec,'all'); % cents
    network_fees = sum((Lambda_b-Lambda_s).*(Es),'all'); % cents
end
end

