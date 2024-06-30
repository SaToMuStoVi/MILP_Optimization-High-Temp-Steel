classdef electrolyzer
    % Class for a electrolyzer unit
    %  properties:
        % OCC;
        % eta;
        % LT;
        % LTP;
        % ACC;
        % fOpex;
        % vOpex;
        % pRatio
        % maxX = 3e4;
        % minX = 0;
        % minLoad = 0.2;                    % min. load in share of capacity
        % maxLoad = 1;                      % max. load in share of capacity
        % Csw

    properties          % parameters
        OCC;                                    % overnight costs €/kW
        eta;                                    % efficiency
        LT;                                     % lifetime in years
        LTP;                                    % lifetime in periods
        ACC;                                    % annualized capital costs €/kW
        fOpex;                                  % cost factor fix OPEX €/kW
        vOpex;                                  % cost factor variable OPEX €/kWh
        pRatio;                                 % pressure ratio
        maxX = 3e4;                             % max. capacity
        minX = 0;                               % min. capacity
        minLoad = 0.2;                          % min. load in share of capacity
        maxLoad = 1;                            % max. load in share of capacity
        Csw;                                    % switching on/off cost factor
        nt;
        ny;
        dt = 1;                                 % 1h timestep
    end
    properties                                  % optimization variables
        Con;                                    % Constraints
        X;                                      % capacity in kW
        acc;
        oper;                                   % component operational variable
        b;                                      % build flag
        r;                                      % retirement flag

        o;                                      % on/off variable
        s;                                      % switching on/off variable
        hOut;                                   % output power in kW
        FC_El;                                  % fuel consumption electricity

        SW_Cost;                                % switching costs
        CAPEX;                                  % Capital investment costs
        FOPEX;                                  % Fix OPEX per period
        VOPEX;                                  % Variable OPEX per period

    end


    methods
        function el = electrolyzer(OCC,eta,LT,fOpex,vOpex,pRatio,maxX,minX,minLoad,maxLoad,Csw,nt,ny)
            % Construct an electrolyzer instance
            %  electrolyzer(OCC,eta,LT,fOpex,vOpex,pRatio,maxX,minX,minLoad,maxLoad,nt,ny)
            if all(OCC>=0) && LT >= 0 && fOpex>=0 && vOpex>=0 && ...
               pRatio>=0 && maxX>=0 && minX>=0 && maxLoad>=0 && minLoad>=0 ...
               && Csw >= 0 
                el.OCC = OCC;
                el.LT = LT;
                el.LTP = LT/5;
                el.ACC = el.OCC/el.LT;
                el.fOpex = fOpex;
                el.vOpex = vOpex;
                el.pRatio = pRatio;
                el.maxX = maxX;
                el.minX = minX;
                el.minLoad = minLoad;
                el.maxLoad = maxLoad;
                el.Csw = Csw;
                el.nt = nt;
                el.ny = ny;
            elseif isa(OCC,'converter')
                el = OCC;
            else
                error('Invalid inputs, all values have to be >= 0');
            end
            if all(0 <= eta)&& all(eta <= 1)
                 el.eta = eta;
            else
                 error('Invalid input: eta either negative or smaller than 1');
            end
            if mod(ny,1) == 0
                el.X = sdpvar(1,ny+el.LTP,'full');                                     % converter capacity in kW
                el.oper = binvar(1,ny+el.LTP,'full');                                  % component operational variable
  
                el.hOut = sdpvar(nt,ny,'full');                                        % output power in kW
                el.b = binvar(1,ny+el.LTP,'full');                                     % build flag
                el.r = binvar(1,ny+el.LTP,'full');                                     % retired flag
                el.o = binvar(nt,ny,'full');                                           % on/off variable
                el.s = sdpvar(nt,ny,'full');                                           % switching on/off variable
                el.FC_El = sdpvar(nt,ny,'full');                                       % electricity consumption in kW
                el.SW_Cost = sdpvar(nt,ny,'full');                                     % switching costs

                el.CAPEX = sdpvar(1,ny,'full');                                        % CAPEX per period
                el.FOPEX = sdpvar(1,ny,'full');                                        % Fix OPEX per period
                el.VOPEX = sdpvar(nt,ny,'full');                                       % Variable OPEX per period
            else
                error('Invalid input, nt has to be a multiple of 5')
            end
            % Design constraints
            el.Con = [el.Con, el.oper*el.minX <= el.X <= el.maxX*el.oper;               % Capacity between min. and max. Capacity, when built 
                
            % Constraints, defining the component lifetime
            el.oper(:,el.LTP+1:end) == el.oper(:,el.LTP:end-1) + el.b(:,el.LTP+1:end)-el.r(:,el.LTP+1:end);   % Component is operational, when built and not retired
            sum(el.b) <= 1;                                                        % Component built just once 
            sum(el.r) <= 1;                                                        % Component retires just once
            el.b(:,1:el.LTP) == 0;                                                 % Component is not built before the modeling horizon
            el.oper(:,1:el.LTP) == 0;                                              % Component is not operational before the modeling horizon
            el.r(:,1:el.LTP) == 0;                                                 % Component is not retired before the modling horizon 
            el.X(:,2:end)-el.X(:,1:end-1) >= 0 + -1e5*(1-el.oper(:,2:end));        % Component has to keep the same capacity during its lifetime
            el.X(:,2:end)-el.X(:,1:end-1) <= 0 + 1e5*(el.b(:,2:end));              % Component has to keep the same capacity during its lifetime
            ];
            
            for i=(el.LTP+1):(el.ny+el.LTP)
            el.Con = [el.Con;
               sum(el.b(:,1:i-(el.LTP))) == sum(el.r(:,1:i));                        % Correlation: component retires, LT years after being built
               ];
            end
            
            % Operational constraints
            el.Con = [el.Con, el.o <= repmat(el.oper(:,el.LTP+1:end),el.nt,1);                     % Component not on, if not operational
                el.o <= repmat(el.oper(:,el.LTP+1:end),nt,1);
                el.hOut >= el.o.*repmat(el.minX,nt,ny);
                el.s(2:el.nt,:) ==  el.o(2:el.nt,:) - el.o(1:el.nt-1,:);                           % switching on/off variable
                el.s(1,:) == el.o(1,:) - el.o(el.nt,:);                                            % switching on/off variable for first and last time step
                el.SW_Cost >= 0; 
                el.SW_Cost >= el.s.*repmat(el.Csw,el.nt,el.ny); 
                el.SW_Cost >= -el.s.*repmat(el.Csw,el.nt,el.ny);                                   % start up costs constraints 
                ];
            
            % Costs
            el.Con = [el.Con, %el.CAPEX == el.ACC.*el.X(:,el.LTP+1:end);                  % Capital costs for one year of every period
                 el.CAPEX >= 0;
                 el.FOPEX == el.fOpex*el.X(:,el.LTP+1:end);                             % Fixed OPEX for one year of every period
                 el.VOPEX == el.vOpex*el.FC_El;                                          % Variable OPEX for the sample week
                ];

            for i=1:el.ny
            el.Con = [el.Con;
                  el.CAPEX >= el.ACC(i)*el.X(:,el.LTP+1:end) - 5e6*(1-el.b(:,i+el.LTP));    % Capital costs: Cost factor fixed to value of the year, the component is built               
               ];
            end

        end
    end
end