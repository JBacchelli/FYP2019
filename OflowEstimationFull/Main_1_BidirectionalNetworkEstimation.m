clc
clear
rng;
%% Initialising fixed variables
[~, A_Bidirectional] = AdjacentMatrix01 (3);%Generate adjacent matrix for n*n network
tau_max = 4;    % maximum path length
n_t = 60;       % observation time horizon
trials = 3;    % number of trials
OD=cell(trials,1);
ODTarget=cell(trials,1);
% Stopping criteria
NMSE_Y_STOP = 5e-4;
for i=1:trials
    %% Data initialize
    [P_initialize,Aod,EList,OList,ODList]=AssignmentMatrix02(A_Bidirectional,tau_max); 
    PC = GetPConstraintsMultiStep01( A_Bidirectional,tau_max,P_initialize );
    [P_target]=AssignmentMatrix02(A_Bidirectional,tau_max);
    n_o=length(OList);
    X_target=OFlowGenerate01(n_t,n_o);
    Y=MVconv02pick(P_target,X_target);
    ODFlow_target=OFlow2ODflow(P_target,X_target,ODList,EList);
    MAX_ITER=250;
    % Initialisation of X (or P), and NMSE_prev
    % P=rand(size(P_initialize));
    P = P_initialize;
    fprintf('P0 has been initalized to target P: %s\n',mat2str(isequal(P_initialize, P_target)));
    NMSE_Y_prev = inf;
        for k=1:1:MAX_ITER
            %% Estimate X
            P_tilde=Toep_P(P,n_t);
            X=Estimate_X(Y,P_tilde,n_o);
            %% Estimate P
            X_tilde=Toep_X(X,tau_max);
            P=Estimate_P0(Y,X_tilde,PC,tau_max);
            %% Calculate ODFlow with estimated X and P
            ODFlow=OFlow2ODflow(P,X,ODList,EList);
            Y_estimated = MVconv02pick(P,X);
            %% Data saving
            ODTarget{i,1}=ODFlow_target;
            NMSE_Y=norm(Y_estimated-Y,'fro')/norm(Y_estimated,'fro');
            history.NMSE_Y(k)=NMSE_Y;
            fprintf('Trial: %d, NMSE_Y: %6.4e, k: %d\n',i,NMSE_Y,k);
            %% Stopping Criteria
            if NMSE_Y < NMSE_Y_STOP
                break
            end
        end
    OD{i,1}=ODFlow;
    save Bi_result OD ODTarget;
end
 %% Dispaly result
[RelativeErrorOD]=TotalODRelativeError(OD,ODTarget);
FigureBar(RelativeErrorOD);