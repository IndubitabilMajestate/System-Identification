%% Figure 1.1
load("product_15.mat");                     % Loading the dataset.
plot(y);                                    % Plotting it for some first impressions.
title('Fig 1.1');
xlabel('t[month]');ylabel('Sales');hold on
%% Figure 1.2
id_size = round(4/5*length(y));             % Splitting the dataset in 80%/20%.
val_size = length(y) - id_size; 
xline(id_size,'k--');
title('Fig 1.2');                           % As stated above we divide the dataset
y_id = y(1:id_size);                        % into 2 parts y_id, y_val for identification
y_val = y(id_size+1:end);                   % and validation respectively. 
t_id = time(1:id_size);                     % Time vector is split as well.
t_val = time(id_size+1:end);  
MSE_id = zeros(1,7);
MSE_val = zeros(1,7);
%% Figure 2.1
figure
plot(y_id);                                 % Plotting y_id for some first impressions.
title('Fig 2.1');
xlabel('t[month]');ylabel('Sales');
%% Figure 2.2(a)(b)
P = 12;                                     % The data is considered to be period with P = 12.
for m = 1:1                                 % The number of Fourier pairs is given by m parameter
    phi_id = zeros(id_size,2*m+2);          % The Phi_id matrix will have the size of id_size*(2*m+2) 
    for i = 1:id_size                       % where m is the current number of Fourier pairs.
        phi_id(i,1) = 1;
        phi_id(i,2) = t_id(i);
        coeff = 1;
        for j = 3:2*m+2
            if mod(j,2) == 1
                phi_id(i,j) = cos(2*coeff*pi*t_id(i)/P);
            else
                phi_id(i,j) = sin(2*coeff*pi*t_id(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    theta = phi_id\y_id;                     % Approximating the theta column vector for m = 1
    y_hat_id = phi_id*theta;                 % We can see that it presents a positive trend along the entire
    MSE_id(m) = mean((y_id-y_hat_id).^2);    % identification dataset, not the best MSE either.
    figure
    plot(y_id);hold; plot(y_hat_id);xlabel('t[month]');ylabel('Sales');                                        
    title('Fig 2.2(a)'); annotation('textbox',[0.2,0.75,0,0],'String',['MSE=',num2str(MSE_id(m))]);

    phi_val = zeros(val_size,2*m+2);         % The Phi_val matrix will have the size of val_size*(2*m+2)
    for i = 1:val_size
        phi_val(i,1) = 1;
        phi_val(i,2) = t_val(i);
        for j = 3:2*m+2
            if mod(j,2) == 1
                phi_val(i,j) = cos(2*coeff*pi*t_val(i)/P);
            else
                phi_val(i,j) = sin(2*coeff*pi*t_val(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    y_hat_val = phi_val*theta;
    MSE_val(m) = mean((y_val-y_hat_val).^2);
    figure
    plot(y_val);hold;plot(y_hat_val);      % Pretty bad results with m = 1 on the validation set,
    title('Fig 2.2(b)');                   % we will calculate the MSE for other values next.
    annotation('textbox',[0.2,0.2,0,0],'String',['MSE=',num2str(MSE_val(m))]);xlabel('t[month]');ylabel('Sales');
end
%% Figure 3
for m = 1:7                                % We will calculate all the MSE_id and MSE_val and plot 
    phi_id = zeros(id_size,2*m+2);         % how the MSE is decreasing/incresing on the validation
    for i = 1:id_size                      % dataset with respect to m.
        phi_id(i,1) = 1;
        phi_id(i,2) = t_id(i);
        coeff = 1;
        for j = 3:2*m+2
            if mod(j,2) == 1
                phi_id(i,j) = cos(2*coeff*pi*t_id(i)/P);
            else
                phi_id(i,j) = sin(2*coeff*pi*t_id(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    theta = phi_id\y_id;                                                    
    y_hat_id = phi_id*theta;                                                
    MSE_id(m) = mean((y_id-y_hat_id).^2);

    phi_val = zeros(val_size,2*m+2);                                        
    for i = 1:val_size
        phi_val(i,1) = 1;
        phi_val(i,2) = t_val(i);
        for j = 3:2*m+2
            if mod(j,2) == 1
                phi_val(i,j) = cos(2*coeff*pi*t_val(i)/P);
            else
                phi_val(i,j) = sin(2*coeff*pi*t_val(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    y_hat_val = phi_val*theta;
    MSE_val(m) = mean((y_val-y_hat_val).^2);
end
figure
plot(MSE_id);hold;plot(MSE_val);
title('Figure 3'); legend('$MSE_{id}$','$MSE_{val}$','Interpreter','latex');xlabel('m');ylabel('MSE');
%% Figures 4.1 & 4.2
m_bestfit = find(MSE_val == min(MSE_val));  % We look for m such that the MSE on the validation 
phi_id = zeros(id_size,2*m_bestfit+2);      % set is minimum and plot the approximation of it.
for i = 1:id_size                                                       
    phi_id(i,1) = 1;
    phi_id(i,2) = t_id(i);
    coeff = 1;
    for j = 3:2*m_bestfit+2
        if mod(j,2) == 1
            phi_id(i,j) = cos(2*coeff*pi*t_id(i)/P);
        else
            phi_id(i,j) = sin(2*coeff*pi*t_id(i)/P);
            coeff = coeff + 1;
        end
    end
end

theta = phi_id\y_id;                                                    
y_hat_id = phi_id*theta;
figure
plot(y_id);hold; plot(y_hat_id);                                        
title('Fig 4.1'); annotation('textbox',[0.2,0.75,0,0],'String',['MSE=',num2str(MSE_id(m_bestfit))]);
xlabel('t[month]');ylabel('Sales');

phi_val = zeros(val_size,2*m_bestfit+2);                                        
for i = 1:val_size
    phi_val(i,1) = 1;
    phi_val(i,2) = t_val(i);
    for j = 3:2*m_bestfit+2
        if mod(j,2) == 1
            phi_val(i,j) = cos(2*coeff*pi*t_val(i)/P);
        else
            phi_val(i,j) = sin(2*coeff*pi*t_val(i)/P);
            coeff = coeff + 1;
        end
    end
end
y_hat_val = phi_val*theta;
figure
plot(y_val);hold;plot(y_hat_val);   % Best result with m = 5 on the validation set% 
title('Fig 4.2');                                                    
annotation('textbox',[0.2,0.2,0,0],'String',['MSE=',num2str(MSE_val(m_bestfit))]);
xlabel('t[month]');ylabel('Sales');
%% Figure 5
y_id = y(1:id_size);
figure
plot(y,'k');hold on                 % Further improvements brought to the reggressors matrix
max_ind = find(y == max(y))+1;      % We find the maximum of the graph and do a piece-wise decomposition
y_id1 = y(1:max_ind);               % Beginning->Maximum (red) has a positive trend
y_id2 = y(max_ind:id_size);         % Maximum->End of identification data (blue) has a negative trend
plot(time(1:max_ind),y_id1,'r',time(max_ind:id_size),y_id2,'b');
xline(max_ind-1,'--k');
title('Figure 5');xlabel('t[month]');ylabel('Sales');
legend('$y_{val}$','$y_{id1}$','$y_{id2}$','$Maximum\, value$','Interpreter','latex','FontSize',15);

%% Figure 6.1 6.2
figure
plot(y_id,'k');
phi_id1 = zeros(max_ind,2);         % We will approximate trends first.
for i = 1:max_ind                   % We use 2 linear regressors for this.
    phi_id1(i,1) = 1;
    phi_id1(i,2) = t_id(i);
end
theta_1 = phi_id1\y_id1;
y_hat_id1 = phi_id1*theta_1;
hold on
plot(y_hat_id1,'b');

phi_id2 = zeros(id_size-max_ind+1,2);
for i = 1:id_size-max_ind+1                                                  
    phi_id2(i,1) = 1;
    phi_id2(i,2) = t_id(i+max_ind-1);
end
theta_2 = phi_id2\y_id2;
y_hat_id2 = phi_id2*theta_2;
diff = y_hat_id1(end)-y_hat_id2(1);
for i = 1:length(y_hat_id2)
    y_hat_id2(i) = y_hat_id2(i) + diff;
end
hold on
plot(t_id(max_ind:id_size),y_hat_id2,'r');
MSE_id(1) = mean((y_id1(1:end-1)-y_hat_id1(1:end-1)).^2) + mean((y_id2-y_hat_id2).^2);
title('Fig 6.1');                                                    
xlabel('t[month]');ylabel('Sales');
for i = 1:max_ind                      % We can subtract the values of each regressor from the data,
    y_id(i) = y_id(i) - y_hat_id1(i);  % thus removing the positive and negative trends.
end
for i = max_ind+1:id_size
    y_id(i) = y_id(i) - y_hat_id2(i-max_ind);
end
figure
plot(y_id);                        % The dataset looks more like some random noise now.
title('Fig 6.2');xlabel('t[month]');ylabel('Sales');legend('$y_{detrend}$','Interpreter','latex');
%% Figure 7
for m = 1:7                        % After detrending the data set we can look for the 
    theta_4 = zeros(2*m+2,1);      % periodicity of the values with the Fourier approximators.
    phi_id3 = zeros(id_size,2*m);
    for i = 1:id_size
        coeff = 1;
        for j = 1:2*m
            if mod(j,2) == 1
                phi_id3(i,j) = cos(2*coeff*pi*t_id(i)/P);
            else
                phi_id3(i,j) = sin(2*coeff*pi*t_id(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    theta_3 = phi_id3\y_id;
    y_hat_id3 = phi_id3*theta_3;
    MSE_id(m) = mean((y_id-y_hat_id3).^2);

    phi_val = zeros(val_size,2*m+2);
    for i = 1:val_size
        coeff = 1;
        phi_val(i,1) = 1;
        phi_val(i,2) = t_val(i);
        for j = 3:2*m+2
            if mod(j,2) == 1
                phi_val(i,j) = cos(2*coeff*pi*t_val(i)/P);
            else
                phi_val(i,j) = sin(2*coeff*pi*t_val(i)/P);
                coeff = coeff + 1;
            end
        end
    end
    theta_4(1) = theta_2(1);
    theta_4(2) = theta_2(2);
    for i = 1:2*m
        theta_4(i+2) = theta_3(i);
    end
    y_hat_val = phi_val*theta_4;
    MSE_val(m) = mean((y_val-y_hat_val).^2);
end
figure
plot(MSE_id);hold;plot(MSE_val);  % The best fit would be m=1, but on the id. dataset its worst.
title('Fig 7');xlabel('m');ylabel('MSE');legend('$MSE_{id}$','$MSE_{val}$','Interpreter','latex');
%% Figure 8
m = 5;                            % We decided to choose m = 5 as the next candidate for the 
theta_4 = zeros(2*m+2,1);         % best approximation and we can see it fits really well on    
phi_id3 = zeros(id_size,2*m);     % the dataset.
for i = 1:id_size
    coeff = 1;
    for j = 1:2*m
        if mod(j,2) == 1
            phi_id3(i,j) = cos(2*coeff*pi*t_id(i)/P);
        else
            phi_id3(i,j) = sin(2*coeff*pi*t_id(i)/P);
            coeff = coeff + 1;
        end
    end
end
theta_3 = phi_id3\y_id;
y_hat_id3 = phi_id3*theta_3;
MSE_id(m) = mean((y_id-y_hat_id3).^2);
phi_val = zeros(val_size,2*m+2);
for i = 1:val_size
    coeff = 1;
    phi_val(i,1) = 1;
    phi_val(i,2) = t_val(i);
    for j = 3:2*m+2
        if mod(j,2) == 1
            phi_val(i,j) = cos(2*coeff*pi*t_val(i)/P);
        else
            phi_val(i,j) = sin(2*coeff*pi*t_val(i)/P);
            coeff = coeff + 1;
        end
    end
end
theta_4(1) = theta_2(1);
theta_4(2) = theta_2(2);
for i = 1:2*m
    theta_4(i+2) = theta_3(i);
end
y_hat_val = phi_val*theta_4;
MSE_val(m) = mean((y_val-y_hat_val).^2);

figure
plot(y_id);hold;plot(y_hat_id3);
title('Fig 8.1');                                                    
annotation('textbox',[0.2,0.2,0,0],'String',['MSE=',num2str(MSE_id(m))]);

figure
plot(y_val);hold;plot(y_hat_val);
title('Fig 8.2');                                                    
annotation('textbox',[0.2,0.2,0,0],'String',['MSE=',num2str(MSE_val(m))]);