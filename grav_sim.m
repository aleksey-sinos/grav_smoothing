function [GA_true, GA_height_m] = grav_sim(STN, SAMPLE_LENGTH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Требуется структура с настройками моделирования гравитационной анамалии
% STN.STN.SG_DGDL_ACT; %действительное значение СКО производной аномалии g по l в 1/c^2
% STN.VEL;         % Скорость объекта вдоль траектории  [м/с]
% STN.STN.SG_GRAV_INSTR;  % СКО шума гравиметра  [м/с^2]
% STN.STN.SG_GRAV_MODEL; % СКО аномалии g [м/с^2]
% STN.SG_h0; %Начальное СКО высоты [м]
% STN.SG_VEL0; %Начальное СКО скорости [м/с]
% 
% %Параметры шума измерений
% STN.SG_ALT; % СКО белого шума измерения высоты  [м]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent stn F G X_DIM W_DIM P C H;
if ~isequal(stn,STN)
    [A,B,C,H] = gravmdl(STN);
    [F,G]=discm(A,B,STN.DT,11,0);
    X_DIM = size(A,1); 
    W_DIM = size(B,2);
    
    if ~isfield(STN, 'P0')
        P=getcov(STN,F,G,C);
    else
        P = STN.P0;
    end
    stn = STN;
end

%выделение памяти
x = zeros(X_DIM,SAMPLE_LENGTH);
GA_true = zeros(1,SAMPLE_LENGTH);
GA_height_m = zeros(1,SAMPLE_LENGTH);

%инициаизация алгоритма
x(:,1)=mvnrnd(zeros(X_DIM,1),P);

for j=1:SAMPLE_LENGTH
    N = randn(W_DIM,1);
    GPS_N = randn();
    x(:,j+1)=F*x(:,j)+G*N;
    GA_true(j) = C*x(:,j);                           %Аномалия
    GA_height_m(j) = H*x(:,j)+STN.SG_ALT*GPS_N;      %"Аномальная высота" с учетом шума гравиметра и шумом СНС
end



    