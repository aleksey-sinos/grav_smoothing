function P = getcov(STN,F,G,C)

%Параметры определения установившегося значения матрицы ковариаций
DK_STAB=1e2;
STAB_THRESHOLD=0.01; 
P0=zeros(size(F,1),size(F,1));
Pg=zeros(3,3);
k=0;

while 1
    k=k+1;
    Pg=F(3:5,3:5)*Pg*F(3:5,3:5)'+G(3:5,2)*G(3:5,2)';
    if abs(sqrt(C(3:5)*Pg*C(3:5)')-STN.SG_GRAV_MODEL)/STN.SG_GRAV_MODEL>STAB_THRESHOLD
        stab_flag=0;
    else
        if ~stab_flag
            stab_flag=1;
            k_stab=k;
        end
        if k-k_stab>=DK_STAB
            break
        end
    end
end % while 1
P0(3:5,3:5)=Pg;

P=P0+blkdiag(diag([STN.SG_h0 STN.SG_VEL0].^2),zeros(3)); %Прямая начальная матрица ковариаций 
end