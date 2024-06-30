classdef storage
    % Class for a storage unit, in the form of a battery or a hydrogen
    % storage tank
    %  properties:
        % OCCPower;
        % OCCCap;
        % LT;
        % LTP;
        % ACCPower;
        % ACCCap;
        % fOpex;
        % vOpex;
        % eta;
        % maxLevel = 100;
        % minLevel = 0;
        % capMax;
        % powerMax;

    properties          % parameters
        OCCPower;                               % overnight costs power
        OCCCap;                                 % overnight costs capacity
        LT;                                     % lifetime in years
        LTP;                                    % lifetime in periods
        ACCPower;                               % annualized capital costs power
        ACCCap;                                 % annualized capital costs capacity
        fOpex;                                  % cost factor fix OPEX
        vOpex = 0;                              % cost factor variable OPEX
        eta;                                    % efficiency
        maxLevel = 100;                         % max. charge level
        minLevel = 0;                           % min. charge level
        capMax;                                 % max. capacity
        powerMax;                               % max. power
        nt;
        ny;
        dt = 1;                                 % 1h timestep
    end
    properties          % optimization variables
        Xcap;                                  % storage capacity in kWh
        Xp;                                    % storage power in kW
        b;                                     % build flag
        r;                                     % retired flag
        oper;                                  % component operational variable
        CH;                                    % storage fill level in kWh
        pChar;                                 % storage energy flow in kW
        pDischar;                              % storage energy flow in kW
        Con;                                   % Constraints

        CAPEX;                                 % overall capital investment costs
        CAPEXp;                                % CAPEX for power
        CAPEXc;                                % CAPEX for capacity
        FOPEX;                                 % Fix OPEX per period
        VOPEX;                                 % Variable OPEX per period
    end


    methods
        function st = storage(OCCPower,OCCCap,LT,fOpex,vOpex,eta,maxLevel,minLevel,capMax,powerMax,nt,ny)
            % Construct a storage instance
            % storage(OCCPower,OCCCap,LT,fOpex,vOpex,eta,maxLevel,minLevel,capMax,powerMax)
            if all(OCCPower>=0) && all(OCCCap>=0) && LT >= 0 && fOpex>=0 && vOpex>=0 && ...
                maxLevel>=0 && minLevel>=0 && capMax>=0 &&  powerMax>=0    
                st.OCCPower = OCCPower;
                st.OCCCap = OCCCap;
                st.LT = LT;
                st.LTP = LT/5;
                st.ACCCap = st.OCCCap/st.LT;
                st.ACCPower = st.OCCPower/st.LT;
                st.fOpex = fOpex;
                st.vOpex = vOpex;
                st.maxLevel = maxLevel;
                st.minLevel = minLevel;
                st.capMax = capMax;
                st.powerMax = powerMax;
                st.nt = nt;
                st.ny = ny;
                st.Xcap = sdpvar(1,ny+st.LTP,'full');
            elseif isa(OCCPower,'storage')
                st = OCCPower;
            else
                error('Invalid inputs, all values have to be >= 0');
            end
            if (0 <= eta)&&(eta <= 1)
                st.eta = eta;
            else
                error('Invalid input: eta either negative or smaller than 1');
            end
            if mod(ny,1) == 0
                st.Xcap = sdpvar(1,ny+st.LTP,'full');
                st.Xp = sdpvar(1,ny+st.LTP,'full');                                    % storage power in kW
                st.b = binvar(1,ny+st.LTP,'full');                                     % build flag
                st.r = binvar(1,ny+st.LTP,'full');                                     % retired flag
                st.oper = binvar(1,ny+st.LTP,'full');                                  % component operational variable
                st.CH = sdpvar(nt,ny,'full');                                   % storage fill level in kWh
                st.pChar = sdpvar(nt,ny,'full');                                % storage energy flow in kW
                st.pDischar = sdpvar(nt,ny,'full');                             % storage energy flow in kW
        
                st.CAPEX = sdpvar(1,ny,'full');                                        % CAPEX per period
                st.CAPEXp = sdpvar(1,ny,'full');                                       % CAPEX power per period
                st.CAPEXc = sdpvar(1,ny,'full');                                       % CAPEX capacity per period
                st.FOPEX = sdpvar(1,ny,'full');                                        % Fix OPEX per period
                st.VOPEX = sdpvar(nt,ny,'full');                                       % Variable OPEX per period            
            
            else
                error('Invalid input, nt has to be a multiple of 5')
            end
                % constraints storage
                st.Con = [st.Con, 0 <= st.Xcap <= st.capMax*st.oper,                   % Capacity between min. and max. Capacity, when built
                          0 <= st.Xp <= st.powerMax*st.oper,                   % power between min. and max. power when built
                st.Xp <= st.Xcap,                                              % power has to be lower than capacity 
            
                % Constraints, defining the component lifetime
                st.oper(:,2:end) == st.oper(:,1:end-1) + st.b(:,2:end)-st.r(:,2:end),  % Component is operational, when built and not retired
                st.oper <= st.Xcap;
                sum(st.b) <= 1,                                                        % Component built just once 
                sum(st.r) <= 1,                                                        % Component retires just once
                st.b(:,1:st.LTP) == 0,                                                 % Component is not built before the modeling horizon
                st.oper(:,1:st.LTP) == 0,                                              % Component is not operational before the modeling horizon
                st.r(:,1:st.LTP) == 0,                                                 % Component is not retired before the modeling horizon
                st.Xcap(:,2:end)-st.Xcap(:,1:end-1) >= 0 + -2e5*(1-st.oper(:,2:end)),  % Component can not decrease its capacity during its lifetime
                st.Xcap(:,2:end)-st.Xcap(:,1:end-1) <= 0 + 2e5*(st.b(:,2:end)),        % Component can not increase its capacity during its lifetime
                st.Xp(:,2:end)-st.Xp(:,1:end-1) >= 0 + -2e5*(1-st.oper(:,2:end)),      % Component can not decrease its power during its lifetime
                st.Xp(:,2:end)-st.Xp(:,1:end-1) <= 0 + 2e5*(st.b(:,2:end)),            % Component can not increase its power during its lifetime
                ];
            
            for i=(st.LTP+1):(st.ny+st.LTP)
                st.Con = [st.Con,
                sum(st.b(:,1:i-(st.LTP))) == sum(st.r(:,1:i))];                     % Correlation: component retires, LT years after being built
            end
                % Operational constraints Battery storage
                st.Con = [st.Con, 0 <= st.CH <= repmat((st.oper(:,st.LTP+1:end)*st.capMax),st.nt,1);                            % fill level must be between 0 and capMax(restriction)
                st.CH <= repmat(st.Xcap(:,st.LTP+1:end),st.nt,1);                                                               % fill level below Xcap(final capacity)
                st.CH >= repmat(st.minLevel*st.Xcap(:,st.LTP+1:end) - 1e6*(1-st.oper(:,st.LTP+1:end)),st.nt,1);                 % min. fill level
                0 <= st.pChar <= repmat(st.powerMax*st.oper(:,st.LTP+1:end),st.nt,1);                                           % flow between +/-powerMax(restriction)
                0 <= st.pChar <= repmat(st.Xp(:,st.LTP+1:end),st.nt,1);                                                         % flow between +/-chosen power Xp(final power)
                0 <= st.pDischar <= repmat(st.powerMax*st.oper(:,st.LTP+1:end),st.nt,1);                                        % flow between +/-powerMax(restriction)
                0 <= st.pDischar <= repmat(st.Xp(:,st.LTP+1:end),st.nt,1);                                                      % flow between +/-chosen power Xp(final power)
                st.CH(2:st.nt,:) == st.CH(1:st.nt-1,:) + st.pChar(1:st.nt-1,:)*st.eta*st.dt - st.pDischar(1:st.nt-1,:)*st.dt;   % charge level in kW
                st.CH(end,:) == 0.5*st.Xcap(:,st.LTP+1:end),                                                                    % Charge level at the end of the sample periods    
                %st.pChar(end,:) == 0,
                %st.pDischar(end,:) == 0,
                ];
             for i = 0:3
                st.Con = [st.Con,                                        
                st.CH(1+(st.nt/4)*i,:) == 0.5*st.Xcap(:,st.LTP+1:end)];                             % Charge level at the beginning and end of every week
             end
                % Costs
                st.Con = [st.Con, st.CAPEX == st.CAPEXp + st.CAPEXc;                                % Capital costs for one year of every period
                st.CAPEX >= 0;
                st.FOPEX == st.fOpex*st.Xcap(:,st.LTP+1:end);                                                       % Fixed OPEX for one year of every period
                st.VOPEX == st.vOpex*st.pChar                                                                       % Variable OPEX for the sample week
                ];

             for i=1:st.ny
                st.Con = [st.Con;
                st.CAPEXp >= st.ACCPower(i)*st.Xp(:,st.LTP+1:end) - 3.2e6*(1-st.b(:,i+st.LTP));    % Capital costs: Cost factor fixed to value of the year, the component is built               
                st.CAPEXc >= st.ACCCap(i)*st.Xcap(:,st.LTP+1:end) - 1.6e5*(1-st.b(:,i+st.LTP));    % Capital costs: Cost factor fixed to value of the year, the component is built 
                ];
            end

        end
    end
end