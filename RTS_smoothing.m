%% ���������
%��������� ������������ �������������:

SAMPLE_LENGTH = 30000;      %������ ����������

%��������� �������������:
STN_SIM.SG_DGDL_ACT=mGAL2MTR(2e-3); % �������������� �������� ��� ����������� �������� g �� l � ����/�
STN_SIM.VEL = 80;         % �������� ������� ����� ����������  [�/�]
STN_SIM.SG_GRAV_INSTR=mGAL2MTR(5);  % ��� ���� ����������  [m���]
STN_SIM.SG_GRAV_MODEL=mGAL2MTR(10); % ��� �������� g [m���]
STN_SIM.SG_h0 = 50e-2;    % ��������� ��� ������ [�]
STN_SIM.SG_VEL0 = 5e-2;   % ��������� ��� �������� [�/�]
STN_SIM.SG_ALT = 0.0158;%5e-2;    % ��� ������ ���� ��������� ������  [�]
STN_SIM.DT = 0.1; %�������� �������������
STN_SIM.A_NOISE = 0;      % ������� �������������� ������
STN_SIM.ALF_ACT = 3*20;   % ��������������  �������� ���������� %�����������
STN_SIM.sh = 60e-2;       % ��� ������������ �������� (��� ������) %�����������

% %��������� ����������:
% STN_EST.SG_DGDL_ACT=mGAL2MTR(2e-3);      % �������������� �������� ��� ����������� �������� g �� l � ����/�
% STN_EST.VEL = 50;                        % �������� ������� ����� ����������  [�/�]
% STN_EST.SG_GRAV_INSTR=mGAL2MTR(50);       % ��� ���� ����������  [m���]
% STN_EST.SG_GRAV_MODEL=mGAL2MTR(10);      % ��� �������� g [m���]
% STN_EST.SG_h0 = 50e-2;                   % ��������� ��� ������ [�]
% STN_EST.SG_VEL0 = 5e-2;                  % ��������� ��� �������� [�/�]
% STN_EST.SG_ALT = 0.0158;%5e-2;                % ��� ������ ���� ��������� ������  [�]
% STN_EST.DT = 0.1; %�������� �������������
% STN_EST.A_NOISE = 0;      % ������� �������������� ������
% STN_EST.ALF_ACT = 3*20;   % ��������������  �������� ���������� %�����������
% STN_EST.sh = 60e-2;       % ��� ������������ �������� (��� ������) %�����������
STN_EST = STN_SIM;

[A,B,C,H] = gravmdl(STN_EST); %������������ ������ �������� � ����� � ����������
X_DIM = size(A,1); Z_DIM = size(B,1); W_DIM = size(B,2); 


%��������� ���������
SERIES_LENGTH=11; % ����� ���� ��� ������������� ������ ��������
[F,G]=discm(A,B,DT,SERIES_LENGTH,0); %���������� ������

X_f = zeros(X_DIM,SAMPLE_LENGTH+1); 
X_f_pr = zeros(X_DIM,SAMPLE_LENGTH+1); 
X_f_cov = zeros(X_DIM,X_DIM,SAMPLE_LENGTH+1);
X_f_cov_pr = zeros(X_DIM,X_DIM,SAMPLE_LENGTH+2);

X_s = zeros(X_DIM,SAMPLE_LENGTH+1); 
X_s_cov = zeros(X_DIM,X_DIM,SAMPLE_LENGTH+1);

GA_f = zeros(SAMPLE_LENGTH,1); 
GA_s = zeros(SAMPLE_LENGTH,1); 

RMSE_s = zeros(1,SAMPLE_LENGTH);
RMSE_f = zeros(1,SAMPLE_LENGTH);

%������������� �������� � ���������
[GA_true, mnt] = grav_sim(STN_SIM,SAMPLE_LENGTH); 

% ��������� �������������� �������������� ������ ��� �������� g
X_f_cov(:,:,1) = getcov(STN_EST,F,G,C);

%������ ������
for k = 2:SAMPLE_LENGTH+1
   X_f_cov_pr(:,:,k) = F*X_f_cov(:,:,k-1)*F'+G*eye(size(G,2))*G';
   K = X_f_cov_pr(:,:,k)*H'*(H*X_f_cov_pr(:,:,k)*H'+STN_EST.SG_ALT^2)^-1;
   X_f_pr(:,k) = F*X_f(:,k-1);
   X_f(:,k) = X_f_pr(:,k)+K*(mnt(k-1)-H*X_f_pr(:,k));
   X_f_cov(:,:,k) = (eye(X_DIM)-K*H)*X_f_cov_pr(:,:,k); 
   RMSE_f(k-1) = sqrt(C*X_f_cov(:,:,k)*C'); %��� ����������
   GA_f(k-1) = C*X_f(:,k); %������ ����������
   
end
X_s(:,SAMPLE_LENGTH+1) = X_f(:,SAMPLE_LENGTH+1);
X_s_cov(:,:,SAMPLE_LENGTH+1) = X_f_cov(:,:,SAMPLE_LENGTH+1);
X_f_cov_pr(:,:,SAMPLE_LENGTH+2) = F*X_f_cov(:,:,SAMPLE_LENGTH+1)*F'+G*eye(size(G,2))*G';

%�������� ������ (�����������)
for k = SAMPLE_LENGTH:-1:1
    K = X_f_cov(:,:,k)*F'*(X_f_cov_pr(:,:,k+1))^-1;
    X_s(:,k) = X_f(:,k)+K*(X_s(:,k+1)-X_f_pr(:,k+1));
    X_s_cov(:,:,k) = X_f_cov(:,:,k)-K*(X_f_cov_pr(:,:,k+1)-X_s_cov(:,:,k+1))*K'; 
    RMSE_s(k) = sqrt(C*X_s_cov(:,:,k)*C'); %��� �����������
    GA_s(k) = C*X_s(:,k); %������ �����������
end

figure(1); clf;   %������ �������
hold on; grid on;
plot((1:SAMPLE_LENGTH)/10000*STN_EST.VEL,MTR2mGAL(GA_true),'k','LineWidth',1);
plot((1:SAMPLE_LENGTH)/10000*STN_EST.VEL,MTR2mGAL(GA_f),'r','LineWidth',1);
plot((1:SAMPLE_LENGTH)/10000*STN_EST.VEL,MTR2mGAL(GA_s),'g','LineWidth',1);
legend('�������� ���','������ ����������','������ �����������','Location','NorthEast')
xlabel('��','fontsize',14)
ylabel('�������� ����','fontsize',14)

figure(2); clf; grid on; hold on; %
plot((1:SAMPLE_LENGTH)/10000*STN_EST.VEL,RMSE_s*mGal,'k','LineWidth',1.5);
plot((1:SAMPLE_LENGTH)/10000*STN_EST.VEL,RMSE_f*mGal,'r--','LineWidth',1.5);
legend('��� �����������','��� ����������','Location','NorthEast')
xlabel('��','fontsize',14)
ylabel('����','fontsize',14)




