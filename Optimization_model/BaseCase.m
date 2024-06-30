% Optimisation model for the identification of decarbonization Pathway for a steel-rolling mill
% Considering the following technologies:
%   - Natural Gas furnace
%   - Induction furnace
%   - Furnace Retrofit for H2
%   - Battery storage
%   - Hydrogen storage tank
%   - Alkaline electrolyzer
%   - PEM electrolyzer
%   - SOEC electrolyzer
%   - Carbon capture and storage

%% data import
%electricityCosts = readmatrix('Input_Data/El_costs_4Periods.csv');
solarIrradiation = readmatrix('Input_Data/Solar_Input_4Periods.csv');
elCosts2030 = readmatrix('Input_Data/Leviathan_2030_4Periods.csv')/1000;
elCosts2040 = readmatrix('Input_Data/Leviathan_2040_4Periods.csv')/1000;
elCosts2050 = readmatrix('Input_Data/Leviathan_2050_4Periods.csv')/1000;

%% time steps
% time vector t in operational layer
T0 = 0;
T = 168*4;
nt = 168*4;             % hourly time steps (4 weeks)
dt = T/nt;   
t = (T0:dt:(T-dt))';    % create time vector

% vector of years y in planning layer
y0 = 0;
Y = 6;                  % 2025 until 2055, 6 periods for installation decisions (2025, 2030, 2035, 2040, 2045, 2050)
ny = 6;
dy = Y/ny;
y = (y0:dy:(Y-dy))';

%Initialize Constraints and Objective function 
Con = [];
obj = 0;

%% General variables
Ir = 0.05;                % Interest rate
cCO2 = [0.0725,0.1013,0.13,0.315,0.5, 0.5];               % Carbon tax in €/kgCO2
gCO2 = [6e7,6e7,6e7,6e7,6e7,0];                           % Goal for emissions of CO2eq in kg per year
PVmax = 5e3;				                              % Max. PV capacity at the location in kW

% Discount factor, considering the different payment points of annual
% payments in a period of 5 years
DFperiod = 0;
for i = 0:5-1                                             % Discount factor for 5 years accumulated
    DFperiod = DFperiod + 1./((1+Ir).^(5*(y+1)-i)); 
end
DFperiod =DFperiod';

CAPEX = sdpvar(1,ny,'full');       % CAPEX per period
FOPEX = sdpvar(1,ny,'full');       % Fix OPEX per period
VOPEX = sdpvar(1,ny,'full');       % Variable OPEX per period
FuelCosts = sdpvar(1,ny,'full');       % Fuel and electricity costs per period
CostsEl = sdpvar(1,ny,'full');       % Electricity costs per period
CostsH = sdpvar(1,ny,'full');       % Hydrogen costs per period
CostsNG = sdpvar(1,ny,'full');       % Natural gas costs per period
CO2Costs = sdpvar(1,ny,'full');        % Carbon costs per period
CO2em = sdpvar(1,ny,'full');           % Carbon emissions per period
SwCosts = sdpvar(1,ny,'full');         % Switching costs per period

% Variables for accumulated Solar PV
PV_CAPEX = sdpvar(1,ny,'full');       % Total CAPEX PV per period
PV_FOPEX = sdpvar(1,ny,'full');       % Total Fix OPEX PV per period
PV_VOPEX = sdpvar(nt,ny,'full');       % Total Variable OPEX PV per period

PV_pOut = sdpvar(nt,ny,'full');        % Accumulated PV output power in kW
PV_CU = sdpvar(nt,ny,'full');          % Accumulated renewable curtailment Solar PV in kW

demand = sdpvar(nt,ny,'full');	       % Heating demand

% Auxillary variables
demand_i = sdpvar(nt,ny,'full');
demand_n = sdpvar(nt,ny,'full');
demand_m = sdpvar(nt,ny,'full');
El_Heat = binvar(1,ny,'full');
Gas_Heat = binvar(1,ny,'full');
GH1_help = binvar(1,ny,'full');
GH2_help = binvar(1,ny,'full');
Mix_Heat = binvar(1,ny,'full');

%% Electricity
% variables
Pgrid = sdpvar(nt,ny,'full');                                           % Grid power demand in kW

% parameters
%EF_el = [0.1575,0,0,0,0,0];          % Emission factor of grid electricity in kg/kWh for carbon neutral el. in 2030
EF_el = [0.131,0.088,0.044,0,0,0];    % Emission factor of grid electricity in kg/kWh for carbon neutral el. in 2040
%EF_el = [0.141,0.112,0.084,0.056,0.028]; % Emission factor of grid electricity in kg/kWh for carbon neutral el. in 2050

% grid contraints
Con = [Con,Pgrid>= 0];

%% Natural gas supply
NG_supply = sdpvar(nt,ny,'full');                             % External supply of natural gas

% parameters
EF_NG = [0.115,0.115,0.115,0.115,0.115,0.115];                % Emission factor of natural gas in kg/kWh

% grid contraints
Con = [Con,NG_supply>= 0];

%% Hydrogen supply
H_supply = sdpvar(nt,ny,'full');                              % External supply of hydrogen

% parameters
EF_H = [0,0,0,0,0,0];                                         % Emission factor of hydrogen in kg/kWh

% grid contraints
Con = [Con,H_supply>= 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation of considered technologies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% inductive heating system
% parameter
Ind = converter(800,0.74,30,30,24.81,0,3e4,0,1,nt,ny);

%% Natural gas furnace
% parameter
NGF = converter(0,0.58,35,15,17.24,0,3e4,0,1,nt,ny);
NGFrb = converter(165,0.58,35,35,17.24,0,3e4,0,1,nt,ny);
Con = [Con,NGFrb.b(:,NGFrb.LTP+1:end)<=1-NGF.oper(:,NGF.LTP+1:end)];
%% Furnace Retrofit for H2
% parameter
H2R = converter(64,0.593,35,35,17.24,0,3e4,0,1,nt,ny);
%% battery storage
%parameters
BS_occp = [240,200,190,175,160,160];
BS_occc = [960,800,760,700,640,640];
BS = storage(BS_occc,BS_occp,15,0,0,0.85,100,0,50000,10000,nt,ny);
BSrb = storage(BS_occc,BS_occp,15,0,0,0.85,100,0,50000,10000,nt,ny);
Con = [Con,BSrb.b(:,BSrb.LTP+1:end)<=1-BS.oper(:,BS.LTP+1:end)];                                % Only rebuilt, if old component retired
%% hydrogen storage
% parameters
HS_occp = [0,0,0,0,0,0];
HS_occc = [33,33,33,33,33,33];
pComp_HS = sdpvar(nt,ny,'full');                                % power for hydrogen compression
HS = storage(HS_occp,HS_occc,20,0.015*15.4,0,1,100,0,200000,10000,nt,ny);
HSrb = storage(HS_occp,HS_occc,20,0.015*15.4,0,1,100,0,200000,10000,nt,ny);
pComp_HSrb = sdpvar(nt,ny,'full');                                % power for hydrogen compression
Con = [Con,HSrb.b(:,HSrb.LTP+1:end)<=sum(HS.r(:,HS.LTP+1:end))];                                % Only rebuilt, if old component retired
%% AEL Electrolyzer
% parameter
AEL_occ = [1000,850,700,550,400,200];
%AEL_occ = [500,400,300,250,220,200];
AEL_eff = [0.65,0.68,0.68,0.68,0.69,0.69];
AEL = electrolyzer(AEL_occ,AEL_eff,30,18.8,0,0.035,3e4,100,0.1,1,1,nt,ny);

%% Hydrogen pipeline
H2_occ = [1000,850,700,550,400,200]*40/30;
%H2_occ = [500,400,300,250,220,200]*40/30;
H2_eff = [1,1,1,1,1,1];
H2p = electrolyzer(H2_occ,H2_eff,40,0,0,0,3e4,100,0,1,0,nt,ny);

%% PEM Electrolyzer
% parameter
% PEM_occ = [1400,1150,900,650,400,200];
% PEM_eff = [0.63,0.63,0.65,0.67,0.68,0.68];
% PEM = electrolyzer(PEM_occ,PEM_eff,30,28,0,0.035,3e4,100,0,1,1,nt,ny);
%% SOEC Electrolyzer
% parameter
% SOEC_occ = [2000,1650,1300,950,600,300];
% SOEC_eff = [0.81,0.83,0.83,0.83,0.83,0.83];
% SOEC = electrolyzer(SOEC_occ,SOEC_eff,10,50,0,0.084,3e4,100,0.2,1,1,nt,ny);
% SOECrb = electrolyzer(SOEC_occ,SOEC_eff,10,50,0,0.084,3e4,0,0.2,1,1,nt,ny);
% Con = [Con,SOECrb.b(:,SOECrb.LTP+1:end)<=sum(SOEC.r(:,SOEC.LTP+1:end))];                                % Only rebuilt, if old component retired
%% Carbon Capture and Storage
% parameter
%CCS = ccs(0.9,25,0.0552,40,NG_supply.*repmat(EF_NG,nt,1),nt,ny);

%% Solar PV
% parameter
PV_occ = [876,300,300,280,280,280];
PV1 = solarPV(1,PV_occ,25,9,0,5e3,solarIrradiation,nt,ny);
%PV2 = solarPV(1,PV_occ,25,9,0,2e3,solarIrradiation,nt,ny);
%PV3 = solarPV(1,PV_occ,25,0,0,1e3,solarIrradiation,nt,ny);
%PV4 = solarPV(1,PV_occ,25,0,0,1e3,solarIrradiation,nt,ny);
%PV5 = solarPV(1,PV_occ,25,0,0,1e3,solarIrradiation,nt,ny);

%% Load demand
demand_Ind = readmatrix("Input_Data/demand_Ind_4Periods.csv");
demand_NGF = readmatrix("Input_Data/demand_NGF_4Periods.csv");
demand_Mix = readmatrix("Input_Data/demand_Mix_4Periods.csv");

%% Fuel prices
% Electricity
c_el = [elCosts2030,elCosts2030,elCosts2040,elCosts2040,elCosts2050,elCosts2050];

% Natural gas price
c_NG_y = [0.0254,0.0217,0.0218,0.0219,0.022,0.022];
c_NG = repmat(c_NG_y,nt,1);

% Hydrogen price
c_H_y = [0.149,0.149,0.1355,0.1222,0.1222,0.1222];                   % hydrogen price for the five periods until 2050, based on pessimistic scenario for import from Spain
c_H = repmat(c_H_y,nt,1);
%% General constraints
Con = [Con, %sum(AEL.b(:,AEL.LTP+1:end) + PEM.b(:,PEM.LTP+1:end) + SOEC.b(:,SOEC.LTP+1:end)) <= 1;   % Only one Electrolyzer technology installed                                      
    %CCS.oper(:,CCS.LTP+1:end) <= GH1_help;             % CCS only if also NGF, NGFrb or H2R and therefore not just Induction heating
    %NGFrb.oper(:,NGFrb.LTP+1:end) + NGF.oper(:,NGF.LTP+1:end) <= 1;     % Either NGF or NGF rebuilt

    El_Heat == 0;                                       % Constraint: Induction furnace is not providing all heat by itself
    
% Defining different heat supply conditions:
     % Gas_Heat = Heat from gas fired furnaces (NG or H2)
     % El_Heat = Heat solo from induction furnace
     % Mix_Heat = Heat from gas fired and induction furnace combined
     
    GH1_help >= NGF.oper(:,NGF.LTP+1:end);
    GH1_help >= NGFrb.oper(:,NGFrb.LTP+1:end);
    GH1_help <= NGF.oper(:,NGF.LTP+1:end) + NGFrb.oper(:,NGFrb.LTP+1:end);
    GH2_help >= H2R.oper(:,H2R.LTP+1:end);
    GH2_help >= GH1_help;
    GH2_help <= H2R.oper(:,H2R.LTP+1:end) + GH1_help;

    Gas_Heat <= GH2_help;
    Gas_Heat <= 1-Ind.oper(:,Ind.LTP+1:end);
    Gas_Heat >= GH2_help + (1-Ind.oper(:,Ind.LTP+1:end))-1;

    El_Heat <= Ind.oper(:,Ind.LTP+1:end);
    El_Heat <= (1-GH2_help); 
    El_Heat >= Ind.oper(:,Ind.LTP+1:end)+(1-GH2_help)-1;

    Gas_Heat + El_Heat + Mix_Heat== 1;

    demand_i == repmat(El_Heat,nt,1).*demand_Ind;
    demand_n == repmat(Gas_Heat,nt,1).*demand_NGF;
    demand_m == repmat(Mix_Heat,nt,1).*demand_Mix;
    demand == demand_i + demand_n + demand_m;
    ];    

%% inductive heating system
% constraints inductive heating
Con = [Con,Ind.Con
         Ind.pOut == Ind.eta*Ind.FC_El;              % Power output in kW
         %Ind.b(:,1:Ind.LTP+3) == 0;                 % Only available from 2040
         Ind.FC_NG == 0;
         Ind.FC_H == 0;
         Ind.pOut <= 0.6*demand;                     % Induction furnace is able to supply only max. 60% of the total heat demand
        ];

%% Natural gas furnace
% constraints Natural gas furnace
Con = [Con, NGF.Con;
    NGF.X(:,NGF.LTP+1) == 16000;                                             % NGF existing at the begininng
    NGF.b(:,NGF.LTP+1) == 1;                                                 % NGF existing at the beginning
    NGF.pOut == NGF.eta*NGF.FC_NG;                                           % power output in kW
    NGF.o <= 1-repmat(H2R.oper(:,H2R.LTP+1:end),nt,1);
    NGF.FC_H == 0;
    NGF.FC_El == 0
    ];

Con = [Con, NGFrb.Con;
    NGFrb.pOut == NGFrb.eta*NGFrb.FC_NG;                                          % power output in kW
    NGFrb.o <= 1-repmat(H2R.oper(:,H2R.LTP+1:end),nt,1);
    NGFrb.FC_H == 0;
    NGFrb.FC_El == 0;
    ];

%% Furnace Retrofit for H2
% constraints H2 Retrofit
Con = [Con, H2R.Con;
     H2R.pOut == H2R.eta*(H2R.FC_NG + H2R.FC_H);              % power output in kW
     H2R.b(:,1:H2R.LTP+1) == 0;                               % Only available from 2030
     H2R.FC_NG >= 0;
     H2R.FC_H >= 0;
     H2R.FC_El == 0;
     H2R.pOut >= 0
     ];

%% battery storage
Con = [Con,BS.Con];
Con = [Con,BSrb.Con];
%% hydrogen storage tank
% constraints hydrogen storage
Con = [Con,HS.Con];
Con = [Con,...
    pComp_HS >= HS.pChar*dt*AEL.pRatio - 1e5*(1-sum(AEL.b));      % power consumption for compression, when AEL
    %pComp_HS >= HS.pChar*dt*PEM.pRatio - 1e5*(1-sum(PEM.b));      % power consumption for compression, when PEM
    %pComp_HS >= HS.pChar*dt*SOEC.pRatio - 1e5*(1-sum(SOEC.b));    % power consumption for compression, when SOEC
    pComp_HS >= 0
    ];

Con = [Con,HSrb.Con];
Con = [Con,...
    pComp_HSrb >= HSrb.pChar*dt*AEL.pRatio - 1e5*(1-sum(AEL.b));      % power consumption for compression, when AEL
    %pComp_HSrb >= HSrb.pChar*dt*PEM.pRatio - 1e5*(1-sum(PEM.b));      % power consumption for compression, when PEM
    %pComp_HSrb >= HSrb.pChar*dt*SOEC.pRatio - 1e5*(1-sum(SOEC.b));    % power consumption for compression, when SOEC
    pComp_HSrb >= 0;
    ];
%% Alkaline Electrolyzer
% Design constraints
Con = [Con,AEL.Con,...
    AEL.o*AEL.minX <= AEL.FC_El <= AEL.maxX*AEL.o;                  % power between 0 and inf when on
    AEL.FC_El <= repmat(AEL.X(:,AEL.LTP+1:end),AEL.nt,1);                       % power below chosen capacity
    AEL.hOut <= repmat(AEL.X(:,AEL.LTP+1:end),AEL.nt,1);
    AEL.hOut == repmat(AEL.eta,nt,1).*AEL.FC_El
    ];
%% Hydrogen pipeline
% Design constraints
Con = [Con,H2p.Con,...
    H2p.o*H2p.minX <= H2p.hOut <= H2p.maxX*H2p.o;                  % power between 0 and inf when on
    H2p.hOut <= repmat(H2p.X(:,H2p.LTP+1:end),H2p.nt,1);           % hydrogen output smaller than capacity
    ];    

%% PEM and SOEC were excluded, as results showed AEL to be the prefered electrolyser for this application

%% PEM Electrolyzer
% Design constraints
% Con = [Con,PEM.Con;
%         PEM.o*PEM.minX <= PEM.FC_El <= PEM.maxX*PEM.o;                  % power between 0 and inf when on
%         PEM.FC_El <= repmat(PEM.X(:,PEM.LTP+1:end),PEM.nt,1);                       % power below chosen capacity
%         PEM.hOut <= repmat(PEM.X(:,PEM.LTP+1:end),PEM.nt,1);
%         PEM.hOut == repmat(PEM.eta,nt,1).*PEM.FC_El
% ];


%% SOEC Electrolyzer
% % Design constraints
% Con = [Con,SOEC.Con;
%         SOEC.o*SOEC.minX <= SOEC.FC_El <= SOEC.maxX*SOEC.o;                  % power between 0 and inf when on
%         SOEC.FC_El <= repmat(SOEC.X(:,SOEC.LTP+1:end),nt,1);            % power below chosen capacity        
%         SOEC.FC_El >= SOEC.minLoad*repmat(SOEC.X(:,SOEC.LTP+1:end),nt,1) - 1e5*repmat((1-SOEC.oper(:,SOEC.LTP+1:end)),nt,1);    % min. load 
%         SOEC.FC_El >= 0;
%         SOEC.hOut == repmat(SOEC.eta,nt,1).*(SOEC.FC_El);                                          % hydrogen output in kW    
%         ];

% Con = [Con,SOECrb.Con;
%         SOECrb.o*SOECrb.minX <= SOECrb.FC_El <= SOECrb.maxX*SOECrb.o;                  % power between 0 and inf when on
%         SOECrb.FC_El <= repmat(SOECrb.X(:,SOECrb.LTP+1:end),nt,1);            % power below chosen capacity        
%         SOECrb.FC_El >= SOECrb.minLoad*repmat(SOECrb.X(:,SOECrb.LTP+1:end),nt,1) - 1e5*repmat((1-SOEC.oper(:,SOEC.LTP+1:end)),nt,1);    % min. load 
%         SOECrb.FC_El >= 0;
%         SOECrb.hOut == repmat(SOECrb.eta,nt,1).*(SOECrb.FC_El);                                          % hydrogen output in kW    
%         ];
%% CCS 
% Design constraints
%Con = [Con,CCS.Con];

%% Solar PV
% Design constraints
Con = [Con,PV1.Con];
%Con = [Con,PV2.Con];
%Con = [Con,PV3.Con];
%Con = [Con,PV4.Con];
%Con = [Con,PV5.Con];
Con = [Con,PV_CAPEX == PV1.CAPEX;% + PV2.CAPEX;% + PV3.CAPEX + PV4.CAPEX + PV5.CAPEX;
          PV_FOPEX == PV1.FOPEX;% + PV2.FOPEX;% + PV3.FOPEX + PV4.FOPEX + PV5.FOPEX;
          PV_VOPEX == PV1.VOPEX;% + PV2.VOPEX;% + PV3.VOPEX + PV4.VOPEX + PV5.VOPEX;
          PV_pOut == PV1.pOut;% + PV2.pOut;% + PV3.pOut + PV4.pOut + PV5.pOut;
          PV_CU == PV1.CU;% + PV2.CU;% + PV3.CU + PV4.CU + PV5.CU
          %repmat(PVmax,1,ny) >= PV1.X(:,PV1.LTP+1:end) + PV2.X(:,PV2.LTP+1:end);                                                                                                                          
   ];

%% Balances
   Con = [Con, Ind.pOut + NGF.pOut + NGFrb.pOut + H2R.pOut == demand;                                   % Heat balance

   Ind.FC_El + AEL.FC_El + PV_CU ...                                                                    % electricity balance
   == Pgrid + BS.pDischar - BS.pChar - pComp_HS + PV_pOut + BSrb.pDischar - BSrb.pChar - pComp_HSrb;

   NGF.FC_NG + NGFrb.FC_NG + H2R.FC_NG == NG_supply;                                                    % Natural gas balance

   H2R.FC_H == H2p.hOut + AEL.hOut + HS.pDischar - HS.pChar ...                                         % Hydrogen balance
   + HSrb.pDischar - HSrb.pChar
   
   H2p.hOut(:,1) == 0;                                                                                  % hydrogen supply only available past 2030
];

%% CO2 emissions constraint
Con = [Con,CO2em == 5*13*sum(Pgrid.*repmat(EF_el,nt,1)+NG_supply.*repmat(EF_NG,nt,1)+H2p.hOut.*repmat(EF_H,nt,1))     % CO2 emissions of operation                                         
        CO2em <= gCO2;                                                                                                 % CO2 emissions smaller than emission goal per period
];

%% Costs
Con = [Con, CAPEX == 5*(Ind.CAPEX + NGF.CAPEX + NGFrb.CAPEX +H2R.CAPEX + ...                          % Capacity Costs
                     BS.CAPEX + HS.CAPEX+ ...
                     AEL.CAPEX + H2p.CAPEX+....
                     PV_CAPEX + BSrb.CAPEX + HSrb.CAPEX);

            FOPEX == 5*(Ind.FOPEX + NGF.FOPEX + NGFrb.FOPEX +H2R.FOPEX + ...                          % Fixed OPEX
                      BS.FOPEX + HS.FOPEX + ...                                          
                      AEL.FOPEX + ...
                      PV_FOPEX + BSrb.FOPEX + HSrb.FOPEX);
    
            VOPEX == 5*13*sum(Ind.VOPEX + NGF.VOPEX + NGFrb.VOPEX + H2R.VOPEX + ...                   % Variable operational costs
                           BS.VOPEX + HS.VOPEX +  ... 
                           AEL.VOPEX + ...
                           PV_VOPEX + BSrb.VOPEX + HSrb.VOPEX);

           CostsEl == 5*13*sum(Pgrid.*c_el);
           CostsH == 5*13*sum(H2p.hOut.*c_H);
           CostsNG == 5*13*sum(NG_supply.*c_NG);
           FuelCosts == CostsEl + CostsH + CostsNG;                                                                                       % Fuel Costs
    
           CO2Costs == 5*13*sum((NG_supply.*repmat(EF_NG,nt,1)+H_supply.*repmat(EF_H,nt,1)).*repmat(cCO2,nt,1)); % -CCS.CO2)+CCS.CO2_Costs);                    % Carbon Costs
    
           SwCosts == 5*13*sum(Ind.SW_Cost + NGF.SW_Cost + NGFrb.SW_Cost + H2R.SW_Cost ...
                          + AEL.SW_Cost);  
           ]

%% objective function
obj = obj + sum(DFperiod.*CAPEX + DFperiod.*FOPEX + DFperiod.*VOPEX + DFperiod.*FuelCosts + DFperiod.*CO2Costs + DFperiod.*SwCosts)/5;


%% optimize
%ops = sdpsettings('verbose',1,'debug',1);
ops = sdpsettings();
ops.verbose = 2;
ops.debug = 2;
ops.showprogress = true;
ops.solver = 'gurobi';
ops.gurobi.MIPGap = 1e-2;
ops.gurobi.TimeLimit = 1e4;
ops.gurobi.SoftMemLimit = 100;

optimize(Con,obj,ops)

%% Output
format short g
costs = value(obj)

% Results Capacities
Cap_NGF = [value(NGF.X(:,NGF.LTP+1:end)),value(NGF.X(:,end))];
Cap_NGFrb = [value(NGFrb.X(:,NGFrb.LTP+1:end)),value(NGFrb.X(:,end))];
Cap_Ind = [value(Ind.X(:,Ind.LTP+1:end)),value(Ind.X(:,end))];
Cap_H2R = [value(H2R.X(:,H2R.LTP+1:end)),value(H2R.X(:,end))];
Cap_BS = [value(BS.Xcap(:,BS.LTP+1:end)),value(BS.Xcap(:,end))];
Cap_BSrb = [value(BSrb.Xcap(:,BSrb.LTP+1:end)),value(BSrb.Xcap(:,end))];
Cap_HS = [value(HS.Xcap(:,HS.LTP+1:end)),value(HS.Xcap(:,end))];
Cap_HSrb = [value(HSrb.Xcap(:,HSrb.LTP+1:end)),value(HSrb.Xcap(:,end))];
Pow_BS = [value(BS.Xp(:,BS.LTP+1:end)),value(BS.Xp(:,end))];
Cap_AEL = [value(AEL.X(:,AEL.LTP+1:end)),value(AEL.X(:,end))];
%Cap_PEM = [value(PEM.X(:,PEM.LTP+1:end)),value(PEM.X(:,end))];
%Cap_SOEC = [value(SOEC.X(:,SOEC.LTP+1:end)),value(SOEC.X(:,end))];
%Cap_CCS = [value(CCS.X(:,CCS.LTP+1:end)),value(CCS.X(:,end))];
 % Cap_PV = [value(PV1.X(:,PV1.LTP+1:end)) + value(PV2.X(:,PV2.LTP+1:end))...
 %     + value(PV3.X(:,PV3.LTP+1:end)) + value(PV4.X(:,PV4.LTP+1:end)) + ...
 %     value(PV5.X(:,PV5.LTP+1:end)), ...
 %     value(PV1.X(:,end)) + value(PV2.X(:,end))... 
 %     + value(PV3.X(:,end)) + value(PV4.X(:,end)) + value(PV5.X(:,end))];
 Cap_PV = [value(PV1.X(:,PV1.LTP+1:end)), value(PV1.X(:,end))];

% Results Operational status
op_NGF = [value(NGF.oper(:,NGF.LTP+1:end)),value(NGF.oper(:,end))];
op_NGFrb = [value(NGFrb.oper(:,NGFrb.LTP+1:end)),value(NGFrb.oper(:,end))];
op_Ind = [value(Ind.oper(:,Ind.LTP+1:end)),value(Ind.oper(:,end))];
op_H2R = [value(H2R.oper(:,H2R.LTP+1:end)),value(H2R.oper(:,end))];
op_BS = [value(BS.oper(:,BS.LTP+1:end)),value(BS.oper(:,end))];
op_HS = [value(HS.oper(:,HS.LTP+1:end)),value(HS.oper(:,end))];
op_BSrb = [value(BSrb.oper(:,BSrb.LTP+1:end)),value(BSrb.oper(:,end))];
op_HSrb = [value(HSrb.oper(:,HSrb.LTP+1:end)),value(HSrb.oper(:,end))];
op_AEL = [value(AEL.oper(:,AEL.LTP+1:end)),value(AEL.oper(:,end))];
%op_PEM = [value(PEM.oper(:,PEM.LTP+1:end)),value(PEM.oper(:,end))];
%op_SOEC = [value(SOEC.oper(:,SOEC.LTP+1:end)),value(SOEC.oper(:,end))];
%op_CCS = [value(CCS.oper(:,CCS.LTP+1:end)),value(CCS.oper(:,end))];
op_PV = [value(PV1.oper(:,PV1.LTP+1:end)),value(PV1.oper(:,end))];

% Costs
capex = value(CAPEX);
fopex = value(FOPEX);
vopex = value(VOPEX);
costsng = value(CostsNG);
costsel = value(CostsEl);
costsh = value(CostsH);
costsco2 = value(CO2Costs);
costsfuel = value(FuelCosts);
costs_el_Ind = 5*13*sum(c_el.*value(Ind.FC_El));
costs_el_AEL = 5*13*sum(c_el.*value(AEL.FC_El));


% Individual costs
indCapex = 5*value(Ind.CAPEX);
ngfCapex = 5*value(NGF.CAPEX);
ngfrbCapex = 5*value(NGFrb.CAPEX);
h2rCapex = 5*value(H2R.CAPEX);
bsCapex = 5*value(BS.CAPEX);
hsCapex = 5*value(HS.CAPEX);
aelCapex = 5*value(AEL.CAPEX);
%pemCapex = 5*value(PEM.CAPEX);
%soecCapex = 5*value(SOEC.CAPEX);
h2pCapex = 5*value(H2p.CAPEX);
pvCapex = 5*value(PV_CAPEX);
bsrbCapex = 5*value(BSrb.CAPEX);
hsrbCapex = 5*value(HSrb.CAPEX);

indFopex = 5*value(Ind.FOPEX);
ngfFopex = 5*value(NGF.FOPEX);
ngfrbFopex = 5*value(NGFrb.FOPEX);
h2rFopex = 5*value(H2R.FOPEX);
bsFopex = 5*value(BS.FOPEX);
hsFopex = 5*value(HS.FOPEX);                                      
aelFopex = 5*value(AEL.FOPEX);
%pemFopex = 5*value(PEM.FOPEX);
%soexFopex = 5*value(SOEC.FOPEX);
pvFopex = 5*value(PV_FOPEX);
bsrbFopex = 5*value(BSrb.FOPEX);
hsrbFopex = 5*value(HSrb.FOPEX);

% Results CO2 Emissions
CO2_El = value(5*13*sum(Pgrid.*repmat(EF_el,nt,1)));
CO2_H = value(5*13*sum(H_supply.*repmat(EF_H,nt,1)));
CO2_NG = value(5*13*sum(NG_supply.*repmat(EF_NG,nt,1)));

% Results Fuel consumption per technology
El_Ind = value(Ind.FC_El);
NG_NGF = value(NGF.FC_NG);
NG_NGFrb = value(NGFrb.FC_NG);
NG_H2R = value(H2R.FC_NG);
H_H2R = value(H2R.FC_H);
El_AEL = value(AEL.FC_El);
H_AEL = value(AEL.hOut);
%El_PEM = value(PEM.FC_El);
%El_SOEC = value(SOEC.FC_El);

% Results storage
pChar_BS = value(BS.pChar);
pDischar_BS = value(BS.pDischar);
pChar_HS = value(HS.pChar);
pDischar_HS = value(HS.pDischar);
pChar_BSrb = value(BSrb.pChar);
pDischar_BSrb = value(BSrb.pDischar);
pChar_HSrb = value(HSrb.pChar);
pDischar_HSrb = value(HSrb.pDischar);

% Results heat output
pOut_Ind = value(Ind.pOut);
pOut_NGF = value(NGF.pOut);
pOut_H2R = value(H2R.pOut);
pOut_NGFrb = value(NGFrb.pOut);

% Fuel consumption
FC_EL_Ind = value(Ind.FC_El);
FC_NG_NGF = value(NGF.FC_NG);
FC_NG_H2R = value(H2R.FC_NG);
FC_NG_NGFrb = value(NGFrb.FC_NG);
FC_H_H2R = value(H2R.FC_H);

% Supply
%H_s = value(H_supply);
P_g = value(Pgrid);
NG_s = value(NG_supply);
H2P_s = value(H2p.hOut);


%% Generate plots
colorNG_NGF = [0.4940, 0.1840, 0.5560]; % purple
colorNG_NGFrb = [0.85, 0.65, 1]; % purple
colorNG_H2R = [0.7, 0.5, 1]; % purple
colorH_H2R = [0.9290, 0.6940, 0.1250];  % yellow
colorH_2 = [1, 0.84, 0];  % yellow
colorEl_Ind = [0, 0.4470, 0.7410];      % blue
colorEl_AEL = [0.53, 0.81, 0.92];      % blue
colorEl_PEM = [0.25, 0.41, 0.88];      % blue
colorEl_SOEC = [0, 0, 0.5];      % blue

color1 = [0.1215, 0.4667,0.7058]; % blau
color2 = [1,0.4980, 0.0549];      % orange
color3 = [0.0902, 0.7451,0.8118]; % türkis
color4 = [0.1725, 0.6274,0.1725]; % grün
color5 = [0.8392, 0.1529,0.1568]; % rot

close all;

% System capacities
figure('Name','Capacities','Position', [100, 100, 800, 600])
p = [2025,2030,2035,2040,2045,2050,2055];
title('Component capacities in kW')
ylabel('capacity in kW');
hold on;
stairs(p,Cap_NGF,'Color',colorNG_NGF);
stairs(p,Cap_NGFrb,'Color',colorNG_NGFrb);
stairs(p,Cap_Ind,'Color',colorEl_Ind);
stairs(p,Cap_H2R,'Color',colorH_H2R);
stairs(p,Pow_BS);
stairs(p,Cap_AEL,'Color',colorEl_AEL);
%stairs(p,Cap_PEM,'Color',colorEl_PEM);
%stairs(p,Cap_SOEC,'Color',colorEl_SOEC);
stairs(p,Cap_PV);
%stairs(p,Cap_CCS);
lgd = legend("Natural gas furnace","NGF rebuilt","Induction furnace", "Hydrogen Retrofit",...
    "Battery storage", "Alkaline Electrolzyer", ...
    "Solar PV",...
    'Location','eastoutside');
hold off

% Storage energy capacities
figure('Name','Storage capacities','Position', [100, 100, 800, 600])
p = [2025,2030,2035,2040,2045,2050,2055];
title('Storage energy capacities in kWh')
ylabel('capacity in kWh');
hold on;
stairs(p,Cap_BS,'color',colorEl_Ind);
stairs(p,Cap_HS,'color',colorH_H2R);
stairs(p,Cap_HSrb,'color',colorH_2);
stairs(p,Cap_BSrb,'color',colorEl_AEL);
lgd = legend("Battery storage", "Hydrogen storage", "HS rebuilt", "BS rebuilt",...
    'Location','eastoutside');
hold off

% Costs
X = {'CAPEX','Fixed OPEX', 'Variable OPEX', 'Costs NG', 'Costs El', 'Costs H' 'Carbon costs'};
bars = [capex; fopex; vopex; costsng; costsel; costsh; costsco2]; 
figure('Name','Costs','Position', [100, 100, 800, 600])
p_costs = ["2025-2030","2030-2035","2035-2040","2040-2045","2045-2050","2050-2055"];
b = bar(p_costs,bars'*10^-6,'stacked');
title('Costs per period')
ylabel('Costs in M€');
lgd = legend(X,'Location', 'northeastoutside');

% CAPEX
X = {'NGF','NGFrb','Ind','H2R','BS','HS','AEL','H2 pipeline','Solar PV','BSrb','HSrb'};
bars = [ngfCapex;ngfrbCapex;indCapex;h2rCapex;bsCapex;hsCapex;aelCapex;h2pCapex;pvCapex;bsrbCapex;hsrbCapex]; 
colors = [colorNG_NGF;colorNG_NGFrb;colorEl_Ind;colorH_H2R;colorEl_PEM;colorH_2;colorEl_AEL;color2;color3;color4;color5];
figure('Name','CAPEX','Position', [100, 100, 800, 600])
p_costs = ["2025-2030","2030-2035","2035-2040","2040-2045","2045-2050","2050-2055"];
b = bar(p_costs,bars'*10^-3,'stacked');
% Set the face color of individual bars
for i = 1:size(bars, 1)
    b(i).FaceColor = colors(i, :);
end
title('Costs per period')
ylabel('Costs in thou.€');
lgd = legend(X,'Location', 'northeastoutside');


% Carbon emissions
X = {'CO2 emissions Electricity','CO2 emissions Natural Gas','CO2 emissions Hydrogen'};
bars = [CO2_El;CO2_NG;CO2_H]; 
colors = [colorEl_Ind;colorNG_NGF;colorH_H2R];
figure('Name','CO2 emissions','Position', [100, 100, 800, 600])
p_costs = ["2025-2030","2030-2035","2035-2040","2040-2045","2045-2050","2050-2055"];
b = bar(p_costs,bars'*10^-6,'stacked');
% Set the face color of individual bars
for i = 1:size(bars, 1)
    b(i).FaceColor = colors(i, :);
end
title('CO2 emissions per period')
ylabel('CO2 emissions in kT');
lgd = legend(X);

% Heat supply by technology
figure('Name','Heat balance','Position', [100, 100, 1200, 800]);
colors = [colorEl_Ind;colorNG_NGF;colorH_H2R;colorNG_NGFrb];
for i=1:6
    subplot(2,3,i)
    data = [pOut_Ind(:,i), value(pOut_NGF(:,i)), value(pOut_H2R(:,i)), value(pOut_NGFrb(:,i))];
    b = bar(data, 'stacked');
    % Set the face color of individual bars
    for j = 1:width(data)
        b(j).FaceColor = colors(j, :);
    end
    hold on 
    plot(c_el(:,i)*1e4,'LineWidth',2,'Color',[0.4660 0.6740 0.1880]);
    plot(H_AEL(:,i),'LineWidth',1,'Color',colorEl_AEL);
    ylim([0, 18000]); % Set Y-axis limits
    xlabel('Time, 4 weeks');
    ylabel('Heat Supply');
end
legend("Induction","Natural gas furnace","Hydrogen Retrofit furnace", "NGF rebuilt","El price in e-4 €","Hydrogen from AEL")

% Fuel and electricity consumption
figure('Name','Fuels','Position', [100, 100, 1000, 800]);
colors = [colorEl_Ind;colorNG_NGF;colorNG_H2R;colorNG_NGFrb;colorH_H2R];
for i=1:6
    subplot(2,3,i)
    data = [value(FC_EL_Ind(:,i)), value(FC_NG_NGF(:,i)), value(FC_NG_H2R(:,i)), value(FC_NG_NGFrb(:,i)),value(FC_H_H2R(:,i))];
    b = bar(data,'stacked');
    % Set the face color of individual bars
    for j = 1:width(data)
        b(j).FaceColor = colors(j, :);
    end
    ylim([0, 18000]); % Set Y-axis limits
    xlabel('Time, 4 weeks');
    ylabel('Fuels');
end
legend("Ind Electricity","NGF natural gas","H2R natural gas", "NGF rebuilt natural gas","H2R hydrogen")

% Fuel consumption per period
X = {'Natural gas, NGF','Natural gas, NGF rebuilt','Natural gas, H2R','Hydrogen, from AEL',...
    'Hydrogen external','Electricity, Ind','Electricity, AEL'};
colors = [colorNG_NGF;colorNG_NGFrb;colorNG_H2R;colorH_H2R;colorH_2;colorEl_Ind;colorEl_AEL];
bars = [5*13*sum(NG_NGF);5*13*sum(NG_NGFrb);5*13*sum(NG_H2R);5*13*sum(H_AEL);5*13*sum(H2P_s);5*13*sum(El_Ind);5*13*sum(El_AEL)]; 
figure('Name','Fuel consumption period','Position', [100, 100, 800, 600])
p_costs = ["2025-2030","2030-2035","2035-2040","2040-2045","2045-2050","2050-2055"];
b = bar(p_costs,bars','stacked');
for i = 1:size(bars, 1)
    b(i).FaceColor = colors(i, :);
end
title('Fuel consumption per period')
ylabel('Fuel consumption in kWh');
lgd = legend(X);

% Component Capacities
figure('Name','Components Capacity','Position', [100, 100, 1000, 800])
p = [2025,2030,2035,2040,2045,2050,2055];
title('Capacity of components')
labels = {'NGF','NGFrb','Ind','H2R','BS','BSrb','HS','HSrb','AEL','PV'};
CapexArray = [Cap_NGF; Cap_NGFrb;Cap_Ind;Cap_H2R;Cap_BS;Cap_BSrb;Cap_HS;Cap_HSrb;Cap_AEL;Cap_PV];
ylimArray = [20000;20000;20000;20000;35000;35000;35000;35000;25000;10000];
for i=1:10
    subplot(10,1,i);
    stairs(p,CapexArray(i,:),'LineWidth', 2);
    ylabel(labels(i));
    ylim([0, ylimArray(i)]); % Set Y-axis limits
end

% Component built binary
figure('Name','Components binary','Position', [100, 100, 1000, 800])
p = [2025,2030,2035,2040,2045,2050,2055];
title('Component built status')
labels = {'NGF','NGFrb','Ind','H2R','BS','BSrb','HS','HSrb','AEL','PV'};
CapexArray = [op_NGF; op_NGFrb;op_Ind;op_H2R;op_BS;op_BSrb;op_HS;Cap_HSrb;op_AEL;op_PV];
for i=1:10
    subplot(10,1,i);
    stairs(p,CapexArray(i,:),'LineWidth', 2);
    ylabel(labels(i));
    ylim([0, 1]); % Set Y-axis limits
end



