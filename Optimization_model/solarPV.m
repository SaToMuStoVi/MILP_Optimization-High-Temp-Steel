classdef solarPV
    % Class for a definite capacity of Solar PV
    %  properties:
        % eta;
        % OCC;
        % LT;
        % LTP;
        % ACC;
        % fOpex;
        % vOpex;
        % Xfix_PV;

    properties
        eta;                % efficiency
        OCC;                % overnight costs in €/kW
        LT;                 % lifetime in years
        LTP;                % lifetime in periods
        ACC;                % annualized capital costs €/kW
        fOpex;              % cost factor fix OPEX €/kW
        vOpex;              % cost factor variable OPEX €/kWh
        Xfix;               % fixed capacity
        solarIrradiation;   % timeseries of solar irradiation during the sample period
        nt;
        ny;
    end
    properties
        Con;                % Constraints of Solar PV
        X;                  % capacity in kW
        oper;               % operational variable
        b;                  % build flag
        r;                  % retirement flag

        pOut;               % output power in kW
        CU;                 % renewable curtailment Solar PV in kW

        CAPEX;              % Capital investment costs
        FOPEX;              % Fix OPEX per period
        VOPEX;              % Variable OPEX per period

    end

    methods
        function pv = solarPV(eta,OCC,LT,fOpex,vOpex,Xfix,solarIrradiation,nt,ny)
            % Construct a solarPV instance
            % solarPV(eta,OCC,LT,fOpex,vOpex,X_fix,nt,ny)
            if all(OCC>=0) && LT >= 0 && fOpex>=0 && vOpex>=0 && Xfix>=0
                pv.OCC = OCC;
                pv.LT = LT;
                pv.LTP = LT/5;
                pv.ACC = pv.OCC/pv.LT;
                pv.fOpex = fOpex;
                pv.vOpex = vOpex;
                pv.Xfix = Xfix;
                pv.solarIrradiation = solarIrradiation;
                pv.nt = nt;
                pv.ny = ny;
            else
                error('Invalid inputs, all values have to be >= 0');
            end 
            if (0 <= eta)&&(eta <= 1)
                pv.eta = eta;
            elseif isa(eta,'solarPV')
                pv = eta;
            else
                error('Invalid input: eta has to be smaller than 1');
            end
            if mod(ny,1) == 0
                pv.pOut = sdpvar(nt,ny,'full');                                        % output power in kW
                pv.CU = sdpvar(nt,ny,'full');                                          % renewable curtailment Solar PV in kW

                pv.X = sdpvar(1,ny+pv.LTP,'full');                                     % converter capacity in kW
                pv.oper = binvar(1,ny+pv.LTP,'full');                                  % component operational variable
                pv.b = binvar(1,ny+pv.LTP,'full');                                     % build flag
                pv.r = binvar(1,ny+pv.LTP,'full');                                     % retired flag

                pv.CAPEX = sdpvar(1,ny,'full');                                        % CAPEX per period
                pv.FOPEX = sdpvar(1,ny,'full');                                        % Fix OPEX per period
                pv.VOPEX = sdpvar(nt,ny,'full');                                       % Variable OPEX per period
            else
                error('Invalid input, nt has to be a multiple of 5')
            end
            % Design constraints
            pv.Con = [pv.Con, pv.X == pv.Xfix*pv.oper;                                 % Capacity fixed, when operational 
                
                % Constraints, defining the component lifetime
                pv.oper(:,pv.LTP+1:end) == pv.oper(:,pv.LTP:end-1) + pv.b(:,pv.LTP+1:end)-pv.r(:,pv.LTP+1:end);   % Component is operational, when built and not retired
                sum(pv.b) <= 1;                                                        % Component built just once 
                sum(pv.r) <= 1;                                                        % Component retires just once
                pv.b(:,1:pv.LTP) == 0;                                                 % Component is not built before the modeling horizon
                pv.oper(:,1:pv.LTP) == 0;                                              % Component is not operational before the modeling horizon
                pv.r(:,1:pv.LTP) == 0;                                                 % Component is not retired before the modling horizon 
                ];
            
            for i=(pv.LTP+1):(pv.ny+pv.LTP)
            pv.Con = [pv.Con;
               sum(pv.b(:,1:i-(pv.LTP))) == sum(pv.r(:,1:i));                          % Correlation: component retires, LT years after being built
               ];
            end
            
            % Operational constraints
            pv.Con = [pv.Con, pv.pOut == pv.eta*(repmat(pv.X(:,pv.LTP+1:end),pv.nt,1).*repmat(pv.solarIrradiation,1,pv.ny)/1000);                      % power output in kW
                    0 <= pv.CU <= pv.pOut;
                ];
            
            % Costs
            pv.Con = [pv.Con, pv.CAPEX >= 0;                  % Capital costs for one year of every period
                 pv.FOPEX == pv.fOpex*pv.X(:,pv.LTP+1:end);                             % Fixed OPEX for one year of every period
                 pv.VOPEX == pv.vOpex*pv.pOut;                                           % Variable OPEX for the sample week
                ];

            for i=1:pv.ny
            pv.Con = [pv.Con;
                  pv.CAPEX >= pv.ACC(i)*pv.X(:,pv.LTP+1:end) - 1.75e5*(1-pv.b(:,i+pv.LTP));    % Capital costs: Cost factor fixed to value of the year, the component is built               
               ];
            end
        end
    end
end