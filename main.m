%--------------------------------------------------------------------------
clc
clear
%--------------------------------------------------------------------------
set(0, 'DefaultAxesFontName', 'Times');
%% grid model

mpc = loadcase('case33re');

%  Active and reactive loads in the retail market: 60% loads in [Baren. Wu.]
mpc.bus(:,3) = mpc.bus(:,3)*1; % active power
mpc.bus(:,4) = mpc.bus(:,4)*0; % reactive power

no = length(mpc.bus);   % # of nodes
br = length(mpc.branch);    % # of branch

fn = mpc.branch(:,1);   % from node
tn = mpc.branch(:,2);   % to node

T = graph(fn,tn);   % graph of network topology
Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);  % make Ybus
Zbus = inv(Ybus(2:end,2:end));    % VSFs
%% network constraints-----change it! true or false

const.activate_Linelimit = true;
const.activate_Voltagelimit = false;
const.activate_Loss = false; % notice : losses are not considered in centralized 

const.Linelimit = [4*ones(1,5) 4*ones(1,6) 1.0*ones(1,21)]';
const.Vmin = 0.95;
const.Vmax = 1.05;

%% Augmented multiplier

rho = 0.02;  % penalty factor

%% For fixed load(nonparticipate p2p energy trade)
%--------------------------------------------------------------------------
Pf = mpc.bus(:,3);  % active power load in the retail market
Qf = mpc.bus(:,4);  % reactive power load in the retail market
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));  % compute power flows
Pline_fix = pfresult.branch(:,14);  % active power flows
Qline_fix = pfresult.branch(:,15);  % reactive power flows
V_fix = pfresult.bus(:,8);  % node voltage

%% Agent environment
[agents, sellers, buyers] = market_info33_new(no, rho);   % market information
% [agents, sellers, buyers] = khorasany(no, rho);   % market information
[centralized_result] = centralized_trading(agents, sellers, buyers, mpc, const);   % centralized_optimization

preference = zeros(length(buyers),length(sellers)); % no preferences
%--------------------------------------------------------------------------
% ADMM trading process (custom setting)
[ADMM.energy, ADMM.sellprice, ADMM.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);  

% first-order gradient method
% [DD.energy, DD.price, DD.LUP] = first_order_trading(sellers, buyers, mpc, const);  

% ADMM trading process (unregulated)
const.activate_Linelimit = false;
const.activate_Voltagelimit = false;
const.activate_Loss = false;
[ADMM_unreg.energy, ADMM_unreg.sellprice, ADMM_unreg.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);

% ADMM trading process (with preference)
preference = 1*rand(length(sellers),length(buyers)); % 0 to 1 cent preference
% preference = zeros(length(buyers),length(sellers)); % no preferences
const.activate_Linelimit = true;
const.activate_Voltagelimit = true;
const.activate_Loss = true;
[ADMM_pref.energy, ADMM_pref.sellprice, ADMM_pref.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);

ISF=makePTDF(mpc);  % injection shift factor (matpower package)
ISF(:,1) = [];  % remove column index corresponding slack bus

% market layer to physical layer mapping
trade_physiclayer = agents.As*sum(ADMM.energy,1)'-agents.Ab*sum(ADMM.energy,2);
trade_physiclayer = trade_physiclayer/1e3; % unit conversion (kW to MW)

trade_physiclayer_unreg = agents.As*sum(ADMM_unreg.energy,1)'-agents.Ab*sum(ADMM_unreg.energy,2);
trade_physiclayer_unreg = trade_physiclayer_unreg/1e3; % unit conversion (kW to MW)

% Assess voltage and line flow
ADMM.VM = [mpc.gen(1,6);V_fix(2:end)+real(Zbus)*(trade_physiclayer)/mpc.baseMVA];
ADMM.Pline = ISF*(trade_physiclayer)+Pline_fix;

ADMM_unreg.VM = [mpc.gen(1,6);V_fix(2:end)+real(Zbus)*(trade_physiclayer_unreg)/mpc.baseMVA];
ADMM_unreg.Pline = ISF*trade_physiclayer_unreg+Pline_fix;

% Assess total trading energy amount
ADMM.sell = sum(ADMM.energy,1)';
ADMM.buy = sum(ADMM.energy,2);

ADMM_unreg.sell = sum(ADMM_unreg.energy,1)';
ADMM_unreg.buy = sum(ADMM_unreg.energy,2);

% DD.sell = sum(DD.energy,2);
% DD.buy = sum(DD.energy,1)';

main_result_graph(T, no, br, agents, sellers, buyers, ADMM, ADMM_unreg, ADMM_pref, centralized_result, mpc, const, preference);

% exactloss = sum(pfresult.branch(:,14)+pfresult.branch(:,16));
% mpc.bus(:,3) = mpc.bus(:,3)-sum(physics_energy,2)/1e3;
% mpc.bus(:,3) = mpc.bus(:,3)+sum(physics_energy,1)'/1e3;
% pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));
% exactloss = sum(pfresult.branch(:,14)+pfresult.branch(:,16));



