function [agent, seller,buyer] = market_info33(no, rho)
% agent description
sellDATA = [
    %  idx bus gen coeffA coeffB genmax
    1   18   0  0.0054  4.34 100
    2   22  40  0.0041  3.27 160
    3   25  10  0.0067  3.69 120
    4   29  10  0.0037  4.24 160
    5   33  10  0.0030  4.35 160];

buyDATA = [
    %  idx bus load coeffA coeffB loadmax
    1   12  30 0.0024  5.49 100
    2   20   10  0.0049  5.21 60
    3   23   10 0.0031  4.89 70
    4   27   40 0.0053  6.04 200
    5   30   40 0.0042  5.85 100];
% a < b/(2*loadmax)
%% struct buyer
buyer(1) = Buyer([1,2,5], buyDATA(1,:), rho);
buyer(2) = Buyer([2,4,5], buyDATA(2,:), rho);
buyer(3) = Buyer([3,4], buyDATA(3,:), rho);
buyer(4) = Buyer([2,3,4], buyDATA(4,:), rho);
buyer(5) = Buyer([1,2,5], buyDATA(5,:), rho);
%% struct seller

for i=1:size(sellDATA,1)
    partner = [];
    for j=1:size(buyDATA,1)
        if ismember(i,buyer(j).partner)
            partner = [partner, j];
        end
    end
    partner = sort(partner);
    seller(i) = Seller(partner, sellDATA(i,:), rho);
end

agent.no = length(seller)+length(buyer);
agent.As = sparse(sellDATA(:,2)',1:length(seller),1,no,length(seller));
agent.Ab = sparse(buyDATA(:,2)',1:length(buyer),1,no,length(buyer));
agent.As(1,:) = []; agent.Ab(1,:) = [];   % remove slack bus
end

