clc
clear
rng;
%% Define folders
STEP = 'single';
COSTS = 'rigid';
ASSMENT = 'random';
H = '3';
W = '3';
TAU_MAX = '4';
TRIALS = '10';
% Folders containing data
folder_data = strcat('./tests/',STEP,'_',COSTS,'_',ASSMENT,'_',H,'by',W,'_',TAU_MAX,'taumax_',TRIALS,'trials/');
% Results folder structured as follows:
%   ./results/date_step_shortpath_sparse_h_w_taumax_nt_trials
SPARSE = 'nosparse';
SAVE = 1;
folder_results = strcat('./results/',STEP,'_',COSTS,'_',ASSMENT,'_',SPARSE,'_',H,'by',W,'_',TAU_MAX,'taumax_',TRIALS,'trials.mat');
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
n_t = 60;       % observation time horizon
n_o = length(o_list);
OD=cell(trials,1);
ODTarget=cell(trials,1);
% Stopping criteria
%DELTA_STOP = 1e-5;
NMSE_Y_STOP = 1e-5;
for i=1:trials
    %% Data initialize
    P_initialise = P_initialise_list(:,:,i);
    P_target = P_target_list(:,:,i);
    %X_target = OFlowGenerate01(n_t,n_o);
    X_target=generate_oflow_single(n_t,n_o);
    Y=P_target * X_target;
    ODFlow_target=oflow2odflow_single(P_target,X_target,od_list,e_list,lc);
    ODTarget{i,1}=ODFlow_target;
    MAX_ITER=1500;
    nmse_history = cell(MAX_ITER,1);
    % Initialisation of X (or P), and NMSE_prev
    % P=rand(size(P_initialize));
    P = P_initialise;
    NMSE_PREV = inf;
    fprintf('P0 has been initalized to target P: %s\n',mat2str(isequal(P_initialise, P_target)));
        for k=1:1:MAX_ITER
            %% Estimate X
            X=estimate_X_single(Y,P,n_t);
            %% Estimate P
            P=estimate_P_single(Y,X,PC,lc);
            %% Calculate ODFlow with estimated X and P
            ODFlow=oflow2odflow_single(P,X,od_list,e_list,lc);
            Y_estimated = P * X;
            NMSE_Y=norm(Y_estimated-Y,'fro')/norm(Y,'fro');%_estimated,'fro');
            DELTA_NMSE = NMSE_PREV - NMSE_Y;
            nmse_history{k,1} = NMSE_Y;
            fprintf('Trial: %d, NMSE_Y: %6.4e, k: %d\n',i,NMSE_Y,k);
            %% Stopping Criteria
            if NMSE_Y < NMSE_Y_STOP %|| DELTA_NMSE < DELTA_STOP
                break
            end
            NMSE_PREV = NMSE_Y;
        end
    OD{i,1}=ODFlow;
end
%% Save all variables of current workspace
% This is needed to keep all information and properly plot data in numpy
if SAVE
    save(folder_results);
end
 %% Dispaly result
[RelativeErrorOD]=TotalODRelativeError(OD,ODTarget);
FigureBarLines(RelativeErrorOD, 1);

function data = load_struct(var, folder)
    data = load(strcat(folder, var, '.mat'));
end

function data = load_var(var, folder)
    load(strcat(folder, var, '.mat'), var);
    data = eval(var);
end