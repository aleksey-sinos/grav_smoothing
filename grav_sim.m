function [GA_true, GA_height_m] = grav_sim(STN, SAMPLE_LENGTH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��������� ��������� � ����������� ������������� �������������� ��������
% STN.STN.SG_DGDL_ACT; %�������������� �������� ��� ����������� �������� g �� l � 1/c^2
% STN.VEL;         % �������� ������� ����� ����������  [�/�]
% STN.STN.SG_GRAV_INSTR;  % ��� ���� ����������  [�/�^2]
% STN.STN.SG_GRAV_MODEL; % ��� �������� g [�/�^2]
% STN.SG_h0; %��������� ��� ������ [�]
% STN.SG_VEL0; %��������� ��� �������� [�/�]
% 
% %��������� ���� ���������
% STN.SG_ALT; % ��� ������ ���� ��������� ������  [�]
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

%��������� ������
x = zeros(X_DIM,SAMPLE_LENGTH);
GA_true = zeros(1,SAMPLE_LENGTH);
GA_height_m = zeros(1,SAMPLE_LENGTH);

%������������ ���������
x(:,1)=mvnrnd(zeros(X_DIM,1),P);

for j=1:SAMPLE_LENGTH
    N = randn(W_DIM,1);
    GPS_N = randn();
    x(:,j+1)=F*x(:,j)+G*N;
    GA_true(j) = C*x(:,j);                           %��������
    GA_height_m(j) = H*x(:,j)+STN.SG_ALT*GPS_N;      %"���������� ������" � ������ ���� ���������� � ����� ���
end



    