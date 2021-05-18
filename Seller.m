classdef Seller
    properties
        bus
        gen
        genmax
        model
        trade_energy
        rho
        partner
        coeff
    end
    
    methods
        function obj = Seller(partner, selldata, rho)
            % initialization
            obj.rho = rho;
            obj.partner = partner;
            obj.bus = selldata(2);
            obj.gen = selldata(3);
            obj.coeff.A = selldata(4);
            obj.coeff.B = selldata(5);
            obj.genmax = selldata(6);
            
            [f, qf, A, b, Aeq, beq, ub, lb] = deal([]);
            
            pos.trade = [0,length(obj.partner)];
            pos.tradesum = [pos.trade(end),pos.trade(end)+1];
            mdl = pos.tradesum(end);       % decision variable for prosumer
            
            tradelb = zeros(1,length(obj.partner));
            tradeub = inf(1,length(obj.partner));
            
            %% Power balance constraint (Equation 1)
            ev = zeros(1,mdl);
            ev(pos.trade(1)+1:pos.trade(end)) = -1; % trade
            ev(pos.tradesum(1)+1:pos.tradesum(end)) = 1; % tradesum
            
            Aeq = [Aeq;ev];
            beq = [beq;0];
            
            % Object function of Prosumer & box constraint
            qf = [qf obj.rho/2*(ones(1,length(obj.partner))) obj.coeff.A];
            f = [f zeros(1,length(obj.partner)) obj.coeff.B];
            lb = [lb tradelb obj.gen];
            ub = [ub tradeub obj.genmax];
            
            obj.model.obj = f;
            obj.model.Q = sparse(diag(qf));
            obj.model.A = sparse([A;Aeq]);
            obj.model.rhs = [b; beq];
            obj.model.sense = [repmat('<',size(A,1),1); repmat('=',size(Aeq,1),1)];
            obj.model.lb = lb;
            obj.model.ub = ub;
        end
        
        function trade_energy = local_opt(obj,eglob,lambda)
            params.outputflag = 0;
            obj.model.obj(1:length(obj.partner)) = obj.rho*(-eglob)-lambda;
            results = gurobi(obj.model, params);
            x = results.x;
            trade_energy = x(1:length(obj.partner));
        end
    end
end

