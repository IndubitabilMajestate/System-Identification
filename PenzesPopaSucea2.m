%% Initialization section
data = load("iddata-20.mat");

%% Testing and plots of the developed model

for m = 1:3
    for na = 1:7
        nb = na;                                               % nb is taken equal to na but testing can be done with na!=nb
        [model1, mse_pred] = nlARX(data.id,data.val,m,na,nb,1);% the mode is set to 1 for prediction results
        [model2, mse_sim] = nlARX(data.id,data.val,m,na,nb,2); % the mode is set to 2 for simulation results
        MSE_id(m,na) = mse_pred(1);
        MSE_val(m,na) = mse_pred(2);

        MSE_sim_val(m,na) = mse_sim(2);
        figure;
        plot(model1.y);hold on; plot(data.val.y);              % we used plot instead of compare due to the instability of
        grid on;title(['MSE=',num2str(mse_pred(2))]);hold off; % the approximation especially on the simulation 
        figure
        plot(model2.y);hold on;plot(data.val.y);
        grid on;title(['MSE=',num2str(mse_sim(2))]);hold off; 
    end
end
%% MSE plots for Simulation and Prediction
for m= 1:3
    figure;
    subplot(1,2,1);plot((1:na),MSE_id(m,:),(1:na),MSE_val(m,:));
    grid on;legend('$MSE_{id}$','$MSE_{val}$','Interpreter','latex','FontSize',15);
    subplot(1,2,2);plot((1:na),MSE_id(m,:),(1:na),MSE_sim_val(m,:));
    grid on;legend('$MSE_{id}$','$MSE_{sim}$','Interpreter','latex','FontSize',15);
end

%% Function Section

function [model,mse] = nlARX(id,val,m,na,nb,mode)           % The function returns a iddata model and the MSE
    u_id = id.u;                                            % for both Identification data and Vaidation data 
    y_id = id.y;
    u_val = val.u;
    y_val = val.y;

    reg_pred_id = RegPredFunc(u_id,y_id,m,na,nb);           % We genrate the regressors with a custom function
    theta = reg_pred_id\y_id;                               % for an easier understanding of the steps
    if mode == 1
        y_pred_id = reg_pred_id*theta;

        mse(1) = mean((y_id-y_pred_id).^2);

        reg_pred_val = RegPredFunc(u_val,y_val,m,na,nb);

        y_pred_val = reg_pred_val*theta;

        mse(2) = mean((y_val-y_pred_val).^2);

        model = iddata(y_pred_val,u_val,val.Ts);            % The prediction model is returned
    elseif mode == 2
        
        y_sim_id = RegSimFunc(u_id,m,na,nb,theta);          % A different function is used for simulation
        y_sim_val = RegSimFunc(u_val,m,na,nb,theta);        % because we generate the output on-the-go
        
        mse(1) = mean((y_id-y_sim_id).^2);
        mse(2) = mean((y_val-y_sim_val).^2);

        model = iddata(y_sim_val,u_val,val.Ts);             % The model is returned at the end(unstable)
    else
        fprintf('Wrong input!\n');                          % Error catching for a different mode
        return
    end
end

function [y_sim] = RegSimFunc(u,m,na,nb,theta)              % Returns the simulated output y
    len = size(powmat(zeros(1,na+nb),m,[]),1);
    y_sim = zeros(length(u),1);
    mat = powmat(zeros(1,na+nb),m,[]);                      % We calculate the power matrix for the elements
    d = zeros(1,na+nb);                                     % in the delay vector after which we will proceed
    for k = 1 : length(u)                                   % to generate the Simulated output one step at a time
        for j = 1:na
            if (k<=j)
                d(j) = 0; 
            else
                d(j) = -y_sim(k-j);
            end
        end
        for j = 1:nb
            if (k<=j)
                d(j+na) = 0; 
            else
                d(j+na) = u(k-j);
            end
        end
        r = genrow(d,mat);                                  % This line generates the regressor row at position k
        y_sim(k) = r*theta;                                 % and then we calculate the output
    end
end

function [Reg] = RegPredFunc(u,y,m,na,nb)                   % Returns the regressor patrix Phi
    len = size(powmat(zeros(1,na+nb),m,[]),1);
    Reg = zeros(length(u),len);
    mat = powmat(zeros(1,na+nb),m,[]);                      % Calculation of the power matrix
    for k = 1 : length(u)
        d = ones(1,na+nb);
        for j = 1:na
            if (k<=j)
                d(j) = 0; 
            else
                d(j) = -y(k-j);
            end
        end
        for j = 1:nb
            if (k<=j)
                d(j+na) = 0; 
            else
                d(j+na) = u(k-j);
            end
        end
        r = genrow(d,mat);                                  % Generating the regressor row at position k
        Reg(k,:) = r;                                       % Adding it to the output matrix
    end
end

function [row] = genrow(d,mat)                              % Returns a vector that contains all the combinations
    row =ones(1,size(mat,1));                               % of the elements of a delay vector risen to their 
    for i = 1:size(mat,1)                                   % respective power from the power matrix
        for j = 1:length(d)
            row(i) = row(i) * d(j)^mat(i,j);
        end
    end
end

function [mat] = powmat(d,m,mat)                            % Returns a matrix that generates all the possible
    if m == 0                                               % power combinations <= to m in a reccursive manner
        mat = zeros(1,length(d));                           % if m=0 returns a vector of 0's which equal to 1
        return                                              % since they represent the power of the elements in d(k)
    else
        mat = powmat(d,m-1,mat);                            % if m>0 we just set the power of one element to m and
        auxd = eye(length(d));                              % subtract from it 1 and distribute it to other elements
        mataux = [];                                        % more like a cascade manner e.g. 
        for i = 1:length(d)                                 % [m 0 0 0]->[m-1 1 0 0]->[m-1 0 1 0]->[m-1 0 0 1]->
            for j = 1:size(mat,1)                           % ->[m-2 2 0 0]->[m-2 1 1 0]->......
                line = mat(j,:) + auxd(i,:);                % Also we check for duplicates, if there are any just to
                f = true;                                   % be sure and add the line to the power matrix
                for k = 1:size(mataux,1)
                    if(mataux(k,:) == line)
                        f = false;
                        break;
                    end
                end
                if(f)
                    mataux = [mataux; line];
                end
            end
        end
        for i = 1:size(mat,1)
            j= 1;
            while j <= size(mataux,1)
                if(mat(i,:) == mataux(j,:))
                    mataux(j,:) = [];
                end
                j = j+1;
            end
        end
        mat = [mat;mataux];
        
    end
end