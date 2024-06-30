classdef converter
    % Class for a converter unit, in the form of a furnace
    %  properties:
        % OCC;
        % eta;
        % LT;
        % RLT;
        % LTP;
        % ACC;
        % fOpex;
        % vOpex;
        % maxX = 3e4;
        % minX = 0;
        % Csw

    properties          % parameters
        OCC;                                    % overnight costs €/kW
        eta;                                    % efficiency/conversion factor
        LT;                                     % lifetime in years
        RLT;                                    % remaining lifetime
        LTP;                                    % lifetime in periods
        ACC;                                    % Annualized capital costs €/kW
        fOpex;                                  % cost factor fix OPEX €/kW
        vOpex;                                  % cost factor variable OPEX €/kWh
        maxX = 3e4;                             % max. capacity
        minX = 0;                               % min capacity
        Csw;                                    % switching on/off cost factor
        nt;
        ny;
        dt = 1;                                 % 1h timestep
    end
    properties                                  % optimization variables
        X;                                      % capacity in kW
        oper;                                   % component operational variable

        pOut;                                   % output power in kW
        b;                                      % build flag
        r;                                      % retirement flag
        o;                                      % on/off variable
        s;                                      % switching on/off variable
        FC_NG;                                  % fuel consumption Natural gas
        FC_H;                                   % fuel consumption hydrogen
        FC_El;                                  % fuel consumption electricity
        Con;                                    % Constraints

        CAPEX;                                  % Capital investment costs
        FOPEX;                                  % Fix OPEX per period
        VOPEX;                                  % Variable OPEX per period
        SW_Cost;                                % switching costs per sample week


    end


    methods
        function cv = converter(OCC,eta,LT,RLT,fOpex,vOpex,maxX,minX,Csw,nt,ny)
            % Construct a converter instance
            % converter(eta,OCC,LT,fOpex,vOpex,maxX,minX,nt,ny)
            if OCC>=0 && LT >= 0 && RLT >=0 && fOpex>=0 && vOpex>=0 && ...
                maxX>=0 && minX>=0 && Csw >=0   
                cv.OCC = OCC;
                cv.LT = LT;
                cv.RLT = RLT;
                cv.LTP = RLT/5;
                cv.ACC = cv.OCC/cv.LT;
                cv.fOpex = fOpex;
                cv.vOpex = vOpex;
                cv.maxX = maxX;
                cv.minX = minX;
                cv.Csw = Csw;
                cv.nt = nt;
                cv.ny = ny;
            elseif isa(OCC,'converter')
                cv = OCC;
            else
                error('Invalid inputs, all values have to be >= 0');
            end
            if (0 <= eta)&&(eta <= 1)
                cv.eta = eta;
            else
                error('Invalid input: eta has to be smaller than 1');
            end
            if mod(ny,1) == 0
                cv.X = sdpvar(1,ny+cv.LTP,'full');                                     % converter capacity in kW
                cv.oper = binvar(1,ny+cv.LTP,'full');                                  % component operational variable

                %cv.pIn = sdpvar(nt,ny,'full');                                        % input power in kW
                cv.pOut = sdpvar(nt,ny,'full');                                        % output power in kW
                cv.b = binvar(1,ny+cv.LTP,'full');                                     % build flag
                cv.r = binvar(1,ny+cv.LTP,'full');                                     % retired flag
                cv.o = binvar(nt,ny,'full');                                           % on/off variable
                cv.s = sdpvar(nt,ny,'full');                                           % switching on/off variable
                cv.FC_NG = sdpvar(nt,ny,'full');                                       % Natural gas consumption in kW
                cv.FC_H = sdpvar(nt,ny,'full');                                        % hydrogen consumption in kW
                cv.FC_El = sdpvar(nt,ny,'full');                                       % electricity consumption in kW
                cv.SW_Cost = sdpvar(nt,ny,'full');                                     % switching costs
                cv.CAPEX = sdpvar(1,ny,'full');                                        % CAPEX per period
                cv.FOPEX = sdpvar(1,ny,'full');                                        % Fix OPEX per period
                cv.VOPEX = sdpvar(nt,ny,'full');                                       % Variable OPEX per period
            else
                error('Invalid input, nt has to be a multiple of 5')
            end
            % Design constraints
            cv.Con = [cv.Con, cv.oper*cv.minX <= cv.X <= cv.maxX*cv.oper;              % Capacity between min. and max. Capacity, when built 
                
                % Constraints, defining the component lifetime
                cv.oper(:,cv.LTP+1:end) == cv.oper(:,cv.LTP:end-1) + cv.b(:,cv.LTP+1:end)-cv.r(:,cv.LTP+1:end);   % Component is operational, when built and not retired
                cv.oper <= cv.X;
                sum(cv.b) <= 1;                                                        % Component built just once 
                sum(cv.r) <= 1;                                                        % Component retires just once
                cv.b(:,1:cv.LTP) == 0;                                                 % Component is not built before the modeling horizon
                cv.oper(:,1:cv.LTP) == 0;                                              % Component is not operational before the modeling horizon
                cv.r(:,1:cv.LTP) == 0;                                                 % Component is not retired before the modeling horizon 
                cv.X(:,2:end)-cv.X(:,1:end-1) >= 0 + -1e5*(1-cv.oper(:,2:end));        % Component has to keep the same capacity during its lifetime
                cv.X(:,2:end)-cv.X(:,1:end-1) <= 0 + 1e5*(cv.b(:,2:end));              % Component has to keep the same capacity during its lifetime
                ];
            
            for i=(cv.LTP+1):(cv.ny+cv.LTP)
            cv.Con = [cv.Con;
               sum(cv.b(:,1:i-(cv.LTP))) == sum(cv.r(:,1:i));                        % Correlation: component retires, LT years after being built
               ];
            end
            
            % Operational constraints
            cv.Con = [cv.Con, cv.o*cv.minX <= cv.pOut <= cv.maxX*cv.o;                  % power between 0 and inf when on
                cv.o <= repmat(cv.oper(:,cv.LTP+1:end),cv.nt,1);                        % Component not on, if not operational
                cv.pOut <= repmat(cv.X(:,cv.LTP+1:end),cv.nt,1);                        % power below chosen capacity
                cv.pOut >= 0;
                cv.s(2:cv.nt,:) ==  cv.o(2:cv.nt,:) - cv.o(1:cv.nt-1,:);                % switching on/off variable
                cv.s(1,:) == cv.o(1,:) - cv.o(cv.nt,:);                                 % switching on/off variable for first and last time step
                cv.SW_Cost >= 0;
                cv.SW_Cost >= cv.s.*repmat(cv.Csw,cv.nt,cv.ny); 
                cv.SW_Cost >= -cv.s.*repmat(cv.Csw,cv.nt,cv.ny);                        % start up costs constraints    
                ];
            
            % Costs
            cv.Con = [cv.Con, cv.CAPEX == cv.ACC*cv.X(:,cv.LTP+1:end);                  % Capital costs for one year of every period
                 cv.FOPEX == cv.fOpex*cv.X(:,cv.LTP+1:end);                             % Fixed OPEX for one year of every period
                 cv.VOPEX == cv.vOpex*cv.pOut;                                           % Variable OPEX for the sample week
                ];
        end
    end
end