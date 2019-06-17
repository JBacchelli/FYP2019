clc
clear
rng;
%% Define folders
STEP = 'multi';
COSTS = 'rigid';
ASSMENT = 'random';
H = '3';
W = '3';
% Folders containing data (add them here as they are created):
%   - tests/multi_any_path
%   - tests/multi_rigid_shortpath
%   - tests/multi_int_shortpath
%   - tests/multi_real_shortpath
folder_data = strcat('./tests/',STEP,'_',COSTS,'_',ASSMENT,'_',H,'by',W,'/');
% Results folder structured as follows:
%   ./results/date_step_shortpath_sparse_h_w_taumax_nt_trials
SAVE = 0;
SPARSE='sparse';
folder_results = strcat('./results/',STEP,'_',COSTS,'_',ASSMENT,'_',SPARSE,'_',H,'by',W,'.mat');
if SAVE && isfile(folder_results)
    disp(strcat('Results in ', folder_results, ' will be overwritten: press any key to confirm'));
    pause;
end

%% Load variables from Python

SHORTEST_PATH = load_var('shortest_path', folder_data);
tau_max = load_var('tau_max', folder_data); % maximum path length
trials = load_var('trials', folder_data);
o_list = load_var('o_list', folder_data);
e_list = load_var('e_list', folder_data);
od_list = load_var('od_list', folder_data);
lc = load_var('lc', folder_data);
P_initialise_list = load_var('P_initialise', folder_data);
P_target_list = load_var('P_target', folder_data);
PC = load_struct('constraints', folder_data);

%% Initialising fixed variables
sparsity_coef = 2; % sparsity coefficient
n_t = 60;       % observation time horizon
n_o = length(o_list);
OD=cell(trials,1);
ODTarget=cell(trials,1);
% Maximum number of iterations
MAX_ITER = 1000;
for i=1:trials
    % Stopping criteria
    L = 3e-4;
    % Data initialize
    P_initialise = P_initialise_list(:,:,:,i);
    P_target = P_target_list(:,:,:,i);
    X_target = OFlowSparseGenerate02(n_o,n_t,sparsity_coef);
    Y = MVconv02pick(P_target, X_target);
    ODFlow_target = oflow2odflow_multi(P_target,X_target,od_list,e_list,lc);
    ODTarget{i,1} = ODFlow_target;
    % Initialisation of X (or P), and NMSE_prev
    % P=rand(size(P_initialize));
    P1 = P_initialise;
    fprintf('P0 has been initalized to target P: %s\n',mat2str(isequal(P_initialise, P_target)));
    %% P1 iteration
    for k = 1:1:MAX_ITER
        % Estimate X
        P_tilde = Toep_P(P1,n_t);
        X = estimate_X_multi(Y,P_tilde,n_o);
        % Estimate P
        X_tilde = Toep_X(X,tau_max);
        P1 = estimate_P_multi(Y,X_tilde,PC,tau_max,n_o,SHORTEST_PATH);
        Y_estimated = MVconv02pick(P1, X);
        NMSE_Y = norm(Y_estimated-Y,'fro')/norm(Y,'fro');%(Y_estimated,'fro');
        history.NMSE_Y(k) = NMSE_Y;
        fprintf('Trial: %d, NMSE_Y: %6.4e, k: %d\n',i,NMSE_Y,k);
        % Stopping Criteria
        if NMSE_Y < L
            break
        end
    end
    %% P2 iteration
    P2 = P1;
    for k2=1:MAX_ITER
        error_bound = L*norm(Y,'fro');
        % Estimate X by l1 norm minimization
        P2_tilde = Toep_P(P2,n_t);
%         X2 = estimate_X_multi_sparsity(Y,P2_tilde,n_o,error_bound);
        X2=estimate_X_multi_sparsity(Y,P2_tilde,n_o,error_bound);
        % Estimate P
        X2_tilde = Toep_X(X2,tau_max);
        P2 = estimate_P_multi(Y,X2_tilde,PC,tau_max,n_o,SHORTEST_PATH);
        % Adjustment process
        ODFlow2 = oflow2odflow_multi(P2,X2,od_list,e_list,lc);
        Y2 = MVconv02pick(P2, X2);
        NMSE_Y2 = norm(Y2-Y,'fro')/norm(Y,'fro')%(Y2,'fro');
        % Calculate sparse transform matrix and L1 norm objective function values
        F = dctmtx(n_t);
        s = repmat('F,',1,n_o);
        Fmtx = eval(sprintf('blkdiag(%s)',s(1:end-1)));
        % Vectorise X
        x2_vec = zeros(1, n_o*n_t);
        for j=1:n_o
            x2_vec(1+n_t*(j-1):j*n_t) = X2(j,:);
        end
        % Value of objective function of L1 minimization at the k2 iteration
        f(k2) = norm(Fmtx*x2_vec',1);
        if k2 >= 2
            fd = f(k2-1)-f(k2);
            fprintf('Trial: %d, k2: %d\n\tfd: %5.3e\n',i,k2,fd);
        end
        % Update ErrorBound
        P2_tilde = Toep_P(P2,n_t);
        X_eb = estimate_X_multi(Y,P2_tilde,n_o);
        Y_eb = MVconv02pick(P2, X_eb);
        curr_eb = norm(Y_eb-Y,'fro');
        fprintf('\tCurrent error: %6.4e, Error bound: %6.4e\n',curr_eb, error_bound);
        if curr_eb <= 0.8*error_bound && k2 > 5 && fd < 1e-2
            L = 0.8*L;
        end
        if error_bound <= 5e-3
            break
        end
     end
    OD{i,1}=ODFlow2;
end
%% Save all variables of current workspace and display results
% This is needed to keep all information and properly plot data in numpy
if SAVE
    save(folder_results);
end
[RelativeErrorOD] = TotalODRelativeError(OD,ODTarget);
FigureBar(RelativeErrorOD, -0.1, 0.1);


function data = load_struct(var, folder)
    data = load(strcat(folder, var, '.mat'));
end

function data = load_var(var, folder)
    load(strcat(folder, var, '.mat'), var);
    data = eval(var);
end