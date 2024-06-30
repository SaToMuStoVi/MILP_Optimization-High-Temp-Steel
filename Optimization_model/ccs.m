classdef ccs
    % Class for a carbon capture and storage unit
    %  properties:
        % lambda;                                 % share saved CO2
        % LT;                                     % lifetime in years
        % LTP;                                    % lifetime in periods
        % cCCS;                                   % Cost per captured CO2: €/tCO2
        % FG;                                     % CO2 amount in Flue gas stream

    properties          % parameters
        lambda;                                 % share saved CO2
        LT;                                     % lifetime in years
        LTP;                                    % lifetime in periods
        cCCS;                                   % Cost per captured CO2: €/tCO2
        cSb;                                    % Standby costs of CCS
        FG;                                     % CO2 amount in Flue gas stream
        nt;
        ny;
        dt = 1;                                 % 1h timestep
    end
    properties                                  % optimization variables
        oper;                                   % component operational variable

        CO2;                                    % saved CO2 in g
        b;                                      % build flag
        r;                                      % retirement flag
        o;                                      % on/off variable
        s;                                      % switching on/off variable
        Con;                                    % Constraints
        CO2_Costs;                              % Costs for Carbon Capture


    end


    methods
        function ccs = ccs(lambda,LT,cCCS,cSb,FG,nt,ny)
            % Construct a CCS instance
            % ccs = ccs(lambda,LT,cCCS,FG,nt,ny)
            if LT >= 0 && cCCS>=0 && cSb>=0 
                ccs.LT = LT;
                ccs.LTP = LT/5;
                ccs.cCCS = cCCS;
                ccs.cSb = cSb;
                ccs.FG = FG;
                ccs.nt = nt;
                ccs.ny = ny;
            elseif isa(lambda,'ccs')
                ccs = lambda;
            else
                error('Invalid inputs, all values have to be >= 0');
            end
            if (0 <= lambda)&&(lambda <= 1)
                ccs.lambda = lambda;
            else
                error('Invalid input: lambda has to be smaller than 1');
            end
            if mod(ny,1) == 0
                ccs.oper = binvar(1,ny+ccs.LTP,'full');                                 % component operational variable

                ccs.CO2 = sdpvar(nt,ny,'full');                                         % saved CO2 in kg
                ccs.b = binvar(1,ny+ccs.LTP,'full');                                    % build flag
                ccs.r = binvar(1,ny+ccs.LTP,'full');                                    % retired flag
                ccs.o = binvar(nt,ny,'full');                                           % on/off variable
                ccs.s = sdpvar(nt,ny,'full');                                           % switching on/off variable
                ccs.CO2_Costs = sdpvar(nt,ny,'full');                                   % Costs for saved CO2
            else
                error('Invalid input, nt has to be a multiple of 5')
            end
            % Design constraints
            ccs.Con = [ccs.Con;
                % Constraints, defining the component lifetime
                ccs.oper(:,2:end) == ccs.oper(:,1:end-1) + ccs.b(:,2:end)-ccs.r(:,2:end); % Component is operational, when built and not retired
                sum(ccs.b) <= 1;                                                          % Component built just once 
                sum(ccs.r) <= 1;                                                          % Component retires just once
                ccs.b(:,1:ccs.LTP) == 0;                                                % Component is not built before the modeling horizon
                ccs.oper(:,1:ccs.LTP) == 0;                                             % Component is not operational before the modeling horizon
                ccs.r(:,1:ccs.LTP) == 0;                                                % Component is not retired before the modeling horizon 
                ];
            
            for i=(ccs.LTP+1):(ccs.ny+ccs.LTP)
            ccs.Con = [ccs.Con;
                sum(ccs.b(:,1:i-(ccs.LTP))) == sum(ccs.r(:,1:i))];                     % Correlation: component retires, LT years after being built
            end
            
            % Operational constraints CCS
            ccs.Con = [ccs.Con, 0 <= ccs.CO2 <= 10e3*ccs.o;                              % CO2 separation between 0 and inf when on
                ccs.o <= repmat(ccs.oper(:,ccs.LTP+1:end),nt,1);
            %0 <= ccs.CO2 <= repmat(10e3.*ccs.oper(:,ccs.LTP+1:end),nt,1);            % CO2 separation between 0 and inf when operational
                ccs.CO2 <= ccs.lambda.*ccs.FG;                                           % CO2 separated output in kgCO2
                ccs.CO2 >= 0;
                ];                                         
            
            % Costs
            ccs.Con = [ccs.Con, ccs.CO2_Costs >= ccs.cCCS.*ccs.CO2                       % Costs of CO2 capture and storage
                ccs.CO2_Costs >= ccs.cSb*repmat(ccs.oper(:,ccs.LTP+1:end),nt,1);
            ];
        end
    end
end