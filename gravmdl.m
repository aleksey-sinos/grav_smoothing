function [A,B,C,H] = gravmdl(STN)

beta=mean(STN.VEL)*(STN.SG_DGDL_ACT)/STN.SG_GRAV_MODEL/sqrt(2);

DZETA=(sqrt(5)-1)/sqrt(5); %Дзета - константа

if STN.A_NOISE == 1
    X_DIM = 6;
else
    X_DIM = 5;
end
W_DIM = 3;
Z_DIM=1; 

A=zeros(X_DIM,X_DIM);
B=zeros(X_DIM,W_DIM);
C=zeros(1,X_DIM);
H=zeros(Z_DIM,X_DIM); H(1,1)=1; 

A(1,2)=1;
A(2,4)=1;
A(3,4)=1;
A(4,5)=1;
A(2,3)=-beta*DZETA;
A(3,3)=-beta;
A(4,4)=-beta;
A(5,5)=-beta;

B(2,1)=STN.SG_GRAV_INSTR;
B(5,2)=STN.SG_GRAV_MODEL*sqrt(10*beta^3);
 
C(1,3)=-beta*DZETA;
C(1,4)=1;

if STN.A_NOISE == 1
    A(6,6) = -1/STN.ALF_ACT;
    B(6,3) = sqrt(2*(1/STN.ALF_ACT)*STN.sh^2);
    H(1,6)=1;
end

%DT = 0.1;
% [F,G]=discm(A,B,DT,11,0);
% 
% % получение установившихся ковариационных матриц для аномалии g
% %Параметры определения установившегося значения матрицы ковариаций
% DK_STAB=1e2;
% STAB_THRESHOLD=0.01;
% P0=zeros(X_DIM,X_DIM);
% Pg=zeros(3,3);
% k=0;
% 
% while 1
%     k=k+1;
%     Pg=F(3:5,3:5)*Pg*F(3:5,3:5)'+G(3:5,2)*G(3:5,2)';
%     if abs(sqrt(C_A(3:5)*Pg*C_A(3:5)')-STN.SG_GRAV_MODEL)/STN.SG_GRAV_MODEL>STAB_THRESHOLD
%         stab_flag=0;
%     else
%         if ~stab_flag
%             stab_flag=1;
%             k_stab=k;
%         end
%         if k-k_stab>=DK_STAB
%             break
%         end
%     end
% end % while 1
% P0(3:5,3:5)=Pg;
% if STN.A_NOISE == 1
%     P=P0+blkdiag(diag([STN.SG_h0 STN.SG_VEL0].^2),zeros(3),STN.sh^2); %Прямая начальная матрица ковариаций всех фильтров банка всех реализаций
% else
%     P=P0+blkdiag(diag([STN.SG_h0 STN.SG_VEL0].^2),zeros(3)); %Прямая начальная матрица ковариаций всех фильтров банка всех реализаций
% end






end