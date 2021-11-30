function [agent, seller,buyer] = market_info33_new(no, rho)
% agent description
sellDATA = [
    %  idx bus gen coeffA coeffB genmax
    1   18  0  0.0058  4.34 220
    2   22  0  0.0038  3.67 260
    3   25  0  0.0027  3.89 180
    4   29  0  0.0034  5.03 240
    5   33  0  0.0040  4.75 160];

buyDATA = [
    %  idx bus load coeffA coeffB loadmax
    1   14  0 0.0024  5.89 100
    2   20   0  0.0042  5.07 180
    3   23   0 0.0031  4.99 180
    4   27   0 0.0021  6.54 260
    5   31   0 0.0018  6.38 240];
% a < b/(2*loadmax)
%% struct buyer
buyer(1) = Buyer([1,2,4], buyDATA(1,:), rho);
buyer(2) = Buyer([1,2,5], buyDATA(2,:), rho);
buyer(3) = Buyer([3,4], buyDATA(3,:), rho);
buyer(4) = Buyer([2,3,4], buyDATA(4,:), rho);
buyer(5) = Buyer([1,3,5], buyDATA(5,:), rho);
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

