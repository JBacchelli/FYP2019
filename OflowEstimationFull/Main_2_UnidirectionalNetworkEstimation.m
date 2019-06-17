clc
clear
rng;
[A_Unidirectional, ~] = AdjacentMatrix01 (3);%[Uni,Bi]
OD=cell(10,1);
ODTarget=cell(10,1);
    for i=1:10%start 100 trials
    %% Data initialize
    tau_max=4;
    tau_max_x=60;
    [P_initialize,Aod,EList,OList,ODList]=AssignmentMatrix02(A_Unidirectional,tau_max); 
    PC = GetPConstraintsMultiStep01( A_Unidirectional,tau_max,P_initialize );
    [P_target]=AssignmentMatrix02(A_Unidirectional,tau_max);
    nn=length(OList);
    X_target=OFlowGenerate01(tau_max_x,nn);
    Y=MVconv02pick(P_target,X_target); 
    ODFlow_target=OFlow2ODflow(P_target,X_target,ODList,EList);
    MAX_ITER=200;
    P=P_initialize;
        for k=1:1:MAX_ITER
        %% Estimate X
        P_tilde=Toep_P(P,tau_max_x);
        X=Estimate_X(Y,P_tilde,nn);
        %% Estimate P
        X_tilde=Toep_X(X,tau_max);
        P=Estimate_P0(Y,X_tilde,PC,tau_max);
        %% Calculate ODFlow with estimated X and P
        ODFlow=OFlow2ODflow(P,X,ODList,EList);
        Y_estimated=MVconv02pick(P,X);
        %% Data saving
        ODTarget{i,1}=ODFlow_target;
        NMSE_Y=norm(Y_estimated-Y,'fro')/norm(Y_estimated,'fro');
        history.NMSE_Y(k)=NMSE_Y;
        fprintf(' NMSE_Y is %e, k is %f \n',NMSE_Y,k);
        %% Stopping Criteria
            if NMSE_Y<8e-4
            break
            end
        end
    OD{i,1}=ODFlow;
    save Uni_result OD ODTarget;
    end
 %% Dispaly result
[RelativeErrorOD]=TotalODRelativeError(OD,ODTarget);
FigureBar(RelativeErrorOD);