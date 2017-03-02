function [F,G]=discm(A,B,dt,N,mode)

n=size(A,1);
r=size(B,2);
Adt=A*dt;
F=eye(n);
for i=N:-1:1
   F=eye(n)+Adt/i*F;
end
if r
   Gm=eye(n);
   for i=N:-1:1
       Gm=eye(n)+Adt*Gm/(i+1);
   end
 
   if mode==0
      G=Gm*B*sqrt(dt);
   else
      G=Gm;
   end
else
   G=[];
end % if r==0   

