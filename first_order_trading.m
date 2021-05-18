function [xij, Lambda, delta_u] = DD_trading_ver2(sellers, buyers, mpc, const)
%% ADMM Energy Trading algorithm with consensus ADMM
%--------------------------------------------------------------------------
% This algorithm use 'Gurobi' optimization solver
% Test feeder data is modified IEEE 34bw
% it was modeled based on matpower case.
%--------------------------------------------------------------------------
%% Network data

no = length(mpc.bus);  % the number of bus
br = length(mpc.branch); % the number of branch

%% Load data read

Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch); % make Ybus
Zbus = inv(Ybus(2:no,2:no))';
ISF=makePTDF(mpc);
%% Energy trading process...

xij = zeros(length(sellers),length(buyers));
yji = zeros(length(sellers),length(buyers));
Lambda = 0*ones(length(sellers),length(buyers));
delta_u = zeros(length(sellers),length(buyers));
delta_l = zeros(length(sellers),length(buyers));

XI_mu = 0.001;
mui_up = zeros(length(sellers),1);
mui_down = zeros(length(sellers),1);
muj_up = zeros(length(buyers),1);
muj_down = zeros(length(buyers),1);
upsilon_up = zeros(br,1);
upsilon_down = zeros(br,1);

MarketPTDF = [];
if const.activate_Linelimit == true
    
    for j=1:length(buyers)
        for i=1:length(sellers)
            MarketPTDF = [MarketPTDF, ISF(:,sellers(i).bus)-ISF(:,buyers(j).bus)];
        end
    end
    MarketPTDF(abs(MarketPTDF)<1e-3)=0;
end

MarketPTDF = MarketPTDF/1e3; % kW to MW

%% For fixed load(nonparticipate p2p energy trade)
%--------------------------------------------------------------------------
Pf = mpc.bus(:,3);
Qf = mpc.bus(:,4);
Pline_fix = -ISF*Pf;
Qline_fix = -ISF*Qf;
V_fix = (-real(Zbus)*Pf(2:end)+imag(Zbus)*Qf(2:end))/mpc.baseMVA;
%--------------------------------------------------------------------------

for round=1:1
    k=1;
    mu_Dual = zeros(br,1);
    ml_Dual = zeros(br,1);
    alpha = 0.1;
    Lambda_prev = 0;
    epsilon_pri = 1;
    epsilon_dual = 1;
    zglob_prev = 0;
    pos_seq = 1;

    while (epsilon_pri > 1e-3) || (epsilon_dual > 1e-2) || (k<10)
        xi_lambda = 0.01/(k^0.009);
        pos_seq = 1/(k^0.01);
        
        for i=1:length(sellers)
            idx_si = sellers(i).partner;
            Lambda(i,idx_si) = max(0,Lambda(i,idx_si)-xi_lambda*(xij(i,idx_si)-yji(i,idx_si))); % eq(21)
            mui_up(i) = max(0,mui_up(i)+XI_mu*(sum(xij(i,idx_si))-sellers(i).genmax)); % eq(22)
            mui_down(i) = max(0,mui_down(i)+XI_mu*(sellers(i).gen-sum(xij(i,idx_si)))); % eq(23)
            xij_hat = (Lambda(i,idx_si)-delta_u(i,idx_si)-mui_up(i)+mui_down(i)-sellers(i).coeff.B)/(2*sellers(i).coeff.A);
            ratiox = (xij(i,idx_si)+pos_seq)/(sum(xij(i,idx_si)+pos_seq));
            xij(i,idx_si) = max(0,xij(i,idx_si)+ratiox.*(xij_hat-sum(xij(i,idx_si)))); % eq(36)
        end
        epsilon_pri = sum(abs(Lambda-Lambda_prev),'all'); % eq(39)
        
        muj_up_prev = muj_up;
        muj_down_prev = muj_down;
        for j=1:length(buyers)
            idx_bj = buyers(j).partner;
            muj_up(j) = max(0,muj_up(j)+XI_mu*(sum(yji(idx_bj,j))-buyers(j).loadmax)); % eq(24)
            muj_down(j) = max(0,muj_down(j)+XI_mu*(buyers(j).load-sum(yji(idx_bj,j)))); % eq(25)
            yji_hat = (-Lambda(idx_bj,j)-delta_u(idx_bj,j)-muj_up(j)+muj_down(j)+buyers(j).coeff.B)/(2*buyers(j).coeff.A);
            ratioy = (yji(idx_bj,j)+pos_seq)/(sum(yji(idx_bj,j)+pos_seq));
            yji(idx_bj,j) = max(0,yji(idx_bj,j)+ratioy.*(yji_hat-sum(yji(idx_bj,j)))); % eq(35)
        end
        epsilon_dual = sum(abs(muj_up-muj_up_prev),'all'); % eq(39)
        epsilon_dual2 = sum(abs(muj_down-muj_down_prev),'all'); % eq(39)
        if const.activate_Linelimit == true
            Pij = MarketPTDF*reshape(xij,[],1)+Pline_fix;
            mu_Dual = max(0,mu_Dual+alpha*(Pij-const.Linelimit));
            ml_Dual = max(0,ml_Dual+alpha*(-const.Linelimit-Pij));
            for i=1:length(sellers)
                for j=1:length(buyers)
                    PTDFij=ISF(:,sellers(i).bus)-ISF(:,buyers(j).bus);
                    delta_u(i,j) = PTDFij'*mu_Dual-PTDFij'*ml_Dual;
                end
            end
        end
        Lambda_prev = Lambda;
        k=k+1;
    end
    utility = 0;
    for j=1:length(buyers)
        utility = utility + (buyers(j).coeff.B*sum(yji(:,j))-buyers(j).coeff.A*sum(yji(:,j))^2);
    end
    for i=1:length(sellers)
        utility = utility - (sellers(i).coeff.B*sum(xij(i,:))+sellers(i).coeff.A*sum(xij(i,:))^2);
    end
    system_fees = sum(delta_u.*(xij),'all');
end
end

