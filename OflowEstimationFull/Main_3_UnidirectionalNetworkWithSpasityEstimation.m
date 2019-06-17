%Main Joint estimation
clc
clear
%close all
rng;
[A_Unidirectional, ~] = AdjacentMatrix01 (3);%[Uni,Bi]
trials = 1;
OD=cell(trials,1);
ODTarget=cell(trials,1);
    %% Starting 10 trials
    for i=1:trials
    %% Data initialization
    f=[];
    tau_max=4;
    [P_initialize,Aod,EList,OList,ODList]=AssignmentMatrix02(A_Unidirectional,tau_max);
    [P_target,~,~,~,~]=AssignmentMatrix02(A_Unidirectional,tau_max);
    PC = GetPConstraintsMultiStep01( A_Unidirectional,tau_max,P_initialize );
    tau_max_x=60;%%%time scale
    nn=length(OList);%%%number of nodes
    S=2;%%%Sparsity density
    [Xtarget,~]=OFlowSparseGenerate02(nn,tau_max_x,S);
    Y=MVconv02pick(P_target,Xtarget);%W_P and W_X help on remove boundary samples 
    ODFlow_target=OFlow2ODflow(P_target,Xtarget,ODList,EList);
    ODTarget{i,1}=ODFlow_target;
    %% Joint Estimation
    MAX_ITER=200;
    P1=P_initialize;
    %% iteration 1 Smart initialization
        for k1=1:1:MAX_ITER
        %% Estimate X
        P_head1=Toep_P(P1,tau_max_x);
        X1=Estimate_X(Y,P_head1,nn);
        %% Estimate P
        X_head1=Toep_X(X1,tau_max);
        P1=Estimate_P0(Y,X_head1,PC,tau_max);
        %% Calculate Y_k1  
        Y_k1=MVconv02pick(P1,X1);
        NMSE_Y1=norm(Y_k1-Y,'fro')/norm(Y_k1,'fro');
        fprintf('Trial: %d, NMSE_Y is %e, k is %f\n',i,NMSE_Y1,k1);
        %%%Stopping criterion
        L=5*1e-4;
            if NMSE_Y1<L
            break 
            end
        end
    %% iteration 2 
    P2=P1;
        for k2=1:1:MAX_ITER
            ErrorBound=L*norm(Y,'fro');
            %% Extimate X by l1 norm minimization
            P_head2=Toep_P(P2,tau_max_x);
            X2=Estimate_X_Sparsity(Y,P_head2,nn,ErrorBound);
            %% Extimate P
            X_head2=Toep_X(X2,tau_max);
            P2=Estimate_P0(Y,X_head2,PC,tau_max);
            %% Adjustment process
            ODFlow2=OFlow2ODflow(P2,X2,ODList,EList);
            Y_k2=MVconv02pick(P2,X2); 
            NMSE_Y2=norm(Y_k2-Y,'fro')/norm(Y,'fro');
            %% Calculate sparse transform matrix and L1 norm objective function values
            F=dctmtx(tau_max_x);
            N = nn;
            s = repmat('F,',1,N);
            Fmtx=eval(sprintf('blkdiag(%s)',s(1:end-1)));%%% Fmtx
            %%% vecterilize X
            X2Vector=[];
                for ii=1:nn
                X2Vector=[X2Vector,X2(ii,:)];
                end
            %%%Value of objective function of L1 minimization at the k2 iteration
            f(k2)=norm(Fmtx*X2Vector',1);
                if k2>=2      
                fd=f(k2-1)-f(k2);
                fprintf('Trial: %d, k2: %d\n\tfd: %5.3e\n',i,k2,fd);
                end
            %% Update ErrorBound
            P_head_eb=Toep_P(P2,tau_max_x);
            X_eb=Estimate_X(Y,P_head_eb,nn);
            Y_eb=MVconv02pick(P2,X_eb);
            curr_eb = norm(Y_eb-Y,'fro');
            fprintf('\tCurrent error: %6.4e, Error bound: %6.4e\n',curr_eb, ErrorBound);
                if curr_eb<=0.5*ErrorBound && k2>5 && fd<1e-2
                        L=0.5*L;       
                end
                if ErrorBound<=1e-6
                break
                end
         end
    OD{i,1}=ODFlow2;
    save SparseUni_result OD ODTarget;
    end
%% Dispaly result
[RelativeErrorOD]=TotalODRelativeError(OD,ODTarget);
FigureBar(RelativeErrorOD);