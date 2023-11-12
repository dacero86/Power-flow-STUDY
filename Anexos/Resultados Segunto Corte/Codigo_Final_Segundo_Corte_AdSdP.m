clc
close all
clear all
%% Toma de datos 

CI = xlsread('Datos.xlsx','Condiciones Iniciales');
Voltajes = CI(:,3);
Angulos = CI(:,4);
Variables =  [Voltajes;Angulos];
EliminacionV = find(isnan(Voltajes(:,1)));
EliminacionO = find(isnan(Angulos(:,1)));
N= length(CI(:,3));
ci = length(EliminacionV);
b = length(EliminacionO);
PPOO = CI(1:3,5:7);
M=0;
for la=1:1:N
    
    do = length(find(isnan(PPOO(la,:))));
    
    if (do==1)
    
        M=M+1;
    else
        
    end    
end


%% codigo de pruebas para otros codigos


% Tramos de linea conectados.
T=[1 3;1 2;2 3];
NL=length(T);

% matruz Ybus
Y=[-i*60 i*40 i*20;i*40 -i*80 i*40;i*20 i*40 -i*60];
Ya=angle(Y);

  
% Numeros de nodos del sistema.
syms V1 V2 V3
syms O1 O2 O3
syms DP1 DP2 DP3
syms DQ1 DQ2 DQ3

V = [V1;V2;V3];
O = [O1;O2;O3];
DP = [DP1; DP2; DP3];
DQ = [DQ1; DQ2; DQ3];

% syms Ya [N 1]


o1=0;
O1=double(o1);
Pesp=[-4 3 -2]';

%% Variables Conocidas 
O1=0;
sol = find(~isnan(Voltajes));
for si=1:1:length(sol)
     
    V(sol(si)) = Voltajes(sol(si));
    
end

fa = find(~isnan(Angulos));

for si=1:1:length(fa)
     
    O(fa(si)) = Angulos(fa(si));
    
end




%% Busqueda de nodos conectados

Busqueda = T(:,1)-T(:,2);

for m = 1:1:N 
    
   
   for k = 1:1:N  
       mi = abs(k-m);
       
       if mi==0
            Pn(k,m) = (abs(V(k))*abs(Y(k,m))*abs(V(m))) * (cos(O(k)-O(m)-Ya(k,m)));
            Qn(k,m) = (abs(V(k))*abs(Y(k,m))*abs(V(m))) * (sin(O(k)-O(m)-Ya(k,m)));
       else
             Conectados = find( (T(:,1)==k & T(:,2)==m)|(T(:,1)==m & T(:,2)==k));
             cu = ~isempty(Conectados);
             if cu==0
                  Pn(k,m) = cu ;
                  Qn(k,m) = cu ;
             else 
                  Pn(k,m) = (abs(V(k))*abs(Y(k,m))*abs(V(m))) * (cos(O(k)-O(m)-Ya(k,m)));
                  Qn(k,m) = (abs(V(k))*abs(Y(k,m))*abs(V(m))) * (sin(O(k)-O(m)-Ya(k,m)));
             end
             
       end 
   end     
end


for aa=1:1:N
    Pkk(aa,1) = sum(Pn(aa,:));   
    Qkk(aa,1) = sum(Qn(aa,:));  
end


%% Eliminacion Incognitas

pot = sym ( 'A' , [b,1]);
Oe = sym ( 'A' , [b,1]);
Pe = sym ( 'A' , [b,1]);
for i=1:1:b
    pot(i) = Pkk(EliminacionO(i));
    Oe(i) = O(EliminacionO(i));
    Pe(i) = DP(EliminacionO(i));
end

poti = sym ( 'A' , [ci,1]);
Ve = sym ( 'A' , [ci,1]);
Qe = sym ( 'A' , [ci,1]);
for i=1:1:ci
    
    poti(ci) = Qkk(EliminacionV(i));
    Ve(ci) = V(EliminacionO(ci));
    Qe(ci) = DQ(EliminacionO(ci));
end

Poteliminado = [pot;poti];
Salidaeliminado = [Oe;Ve];
Entradaelminado = [Pe;Qe];
Ni = length(Salidaeliminado)


%% Jacobiano
Jacobiano = jacobian(Poteliminado,Salidaeliminado);

%% Componentes del Jacobiano 
Hm = Jacobiano(1:(N-1),1:(N-1));
Nm = Jacobiano(1:(N-1),N:(N-1+M));
Jm = Jacobiano(N:(N-1+M),1:(N-1));
Lm = Jacobiano((N-1+M),(N-1+M));
Jacobianom = [Hm,Nm;Jm,Lm];




%% Condiciones iniciales 

re_s = length(EliminacionO);
re_b = length(EliminacionV);


for viv=1:1:b
    
    O(EliminacionO(viv)) = 0;
    
end


viv;

for vivi=1:1:ci
    
    V(EliminacionV(vivi)) = 1;
    
end

 
for k=1:N
eval(sprintf(' O%d = O(k) ', k));
end


for k=1:N
eval(sprintf(' V%d = V(k) ', k));
end

%% Switch case 

disp('Ingrese las iniciales del metodo que desea que se desarrolle')
disp('NR si desea el metodo Newthon Raphson')
disp('NRM si desea el metodo Newthon Raphson Modificado')
disp('NRD si desea el metodo Newthon Raphson Desacoplado')
disp('NRDR si desea el metodo Newthon Raphson Desacoplado Rapido')
METODO = input('Ingrese las iniciales: ');

switch METODO
    case 'NR'
         %% Newthon Raphson
         tic
         nrp = zeros(1,length(Poteliminado));
         nrs = zeros(1,length(Salidaeliminado));
         
         Pcal0=eval(Poteliminado);
         Jacobiano0=eval(Jacobiano);
         DeltaP0=Pesp-Pcal0;
         DeltaP0 = vpa(DeltaP0,6);
         Delta0=(inv(Jacobiano0))*DeltaP0;
         Variable0= eval(Salidaeliminado) + Delta0;
         
         nrp(1,1:length(Poteliminado)) =  DeltaP0;
         nrs(1,1:length(Salidaeliminado)) = Variable0;
         
         e=0.001;
         a=600; 
         L=0; % Marca la iteracion 0
         
         
         for k=1:1:length(EliminacionO)
         eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
         end
         
         
         for k=1:1:length(EliminacionV)
         eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
         end
         
         
         
         
         detecta = 0;
         
         while Ni ~= detecta
              detecta = 0;
              for f = 1:1:Ni
                  t = e <= abs(DeltaP0(f,1));
                  t = eval(t);
                  detecta = ~t + detecta ;
              end
         
              L = L + 1 ;
              Pcal0=eval(Poteliminado);
              Jacobiano0=eval(Jacobiano);
              DeltaP0=Pesp-Pcal0;
              DeltaP0 = vpa(DeltaP0,6);
              Delta0=(inv(Jacobiano0))*DeltaP0;
              Variable0= eval(Salidaeliminado)+(Delta0);
              Variable0 = vpa(Variable0,6);
              nrp(L+1,1:length(Poteliminado)) =  DeltaP0;
              nrs(L+1,1:length(Salidaeliminado)) = Variable0;
            
              for k=1:1:length(EliminacionO)
              eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
              end
         
         
              for k=1:1:length(EliminacionV)
              eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
              end
         
              
         end
         clc
         toc
         disp('Numero de iteraciones:')
         disp(L)
         disp('Delta de potencia por cada iteración:')
         disp([Entradaelminado.';vpa(nrp)])
         disp('Variables Objetivo por cada iteración:')
         disp([Salidaeliminado.';vpa(nrs)])
         
    case 'NRM'
        tic
        nrmp = zeros(1,length(Poteliminado));
        nrms = zeros(1,length(Salidaeliminado));
        
        Nm = Nm*V2;
        Lm = Lm*V2;
        Jacobianom = [Hm,Nm;Jm,Lm];
        Jacobiano = Jacobianom; 
        
        Pcal0=eval(Poteliminado);
        Jacobiano0=eval(Jacobiano);
        DeltaP0=Pesp-Pcal0;
        Delta0=(inv(Jacobiano0))*DeltaP0;
        NRM = length(pot)+length(poti);
        Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
        Variable0= eval(Salidaeliminado)+ Delta0;
        
        nrmp(1,1:length(Poteliminado)) =  DeltaP0;
        nrms(1,1:length(Salidaeliminado)) = Variable0;
        
        e=0.001;
        a=600; 
        L=0; % Marca la iteracion 0
        
        
        for k=1:1:length(EliminacionO)
        eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
        end
        
        
        for k=1:1:length(EliminacionV)
        eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
        end
        
        detecta = 0;
        
        while Ni ~= detecta
             detecta = 0;
             for f = 1:1:Ni
                 t = e <= abs(DeltaP0(f,1));
                 t = eval(t);
                 detecta = ~t + detecta ;
             end
        
             L = L + 1 ;
             Pcal0=eval(Poteliminado);
             Jacobiano0=eval(Jacobiano);
             DeltaP0=Pesp-Pcal0;
             DeltaP0 = vpa(DeltaP0,6);
             Delta0=(inv(Jacobiano0))*DeltaP0;
             NRM = length(pot)+length(poti);
             Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
             Variable0= eval(Salidaeliminado)+(Delta0);
             Variable0 = vpa(Variable0,6);
             nrmp(L+1,1:length(Poteliminado)) =  DeltaP0;
             nrms(L+1,1:length(Salidaeliminado)) = Variable0;
           
             for k=1:1:length(EliminacionO)
             eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
             end
        
        
             for k=1:1:length(EliminacionV)
             eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
             end
        
             
        end
        clc
        toc
        disp('Numero de iteraciones:')
        disp(L)
        disp('Delta de potencia por cada iteración:')
        disp([Entradaelminado.';vpa(nrmp)])
        disp('Variables Objetivo por cada iteración:')
        disp([Salidaeliminado.';vpa(nrms)])
    case 'NRD'
            %% Newthon Raphson desacoplado
            tic
            nrdp = zeros(1,length(Poteliminado));
            nrds = zeros(1,length(Salidaeliminado));

            Nm = Nm*V2;
            Lm = Lm*V2;
            Jm = zeros(size(Jm));
            Nm = zeros(size(Nm));
            Jacobianom = [Hm,Nm;Jm,Lm];
            Jacobiano = Jacobianom;

            Pcal0=eval(Poteliminado);
            Jacobiano0=eval(Jacobiano);
            DeltaP0=Pesp-Pcal0;
            Delta0=(inv(Jacobiano0))*DeltaP0;
            NRM = length(pot)+length(poti);
            Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
            Variable0= eval(Salidaeliminado)+ Delta0;

            nrdp(1,1:length(Poteliminado)) =  DeltaP0;
            nrds(1,1:length(Salidaeliminado)) = Variable0;

            e=0.001;
            a=600; 
            L=0; % Marca la iteracion 0

            for k=1:1:length(EliminacionO)
            eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
            end


            for k=1:1:length(EliminacionV)
            eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
            end

            detecta = 0;

            while Ni ~= detecta

                 detecta = 0;
                 for f = 1:1:Ni
                     t = e <= abs(DeltaP0(f,1));
                     t = eval(t);
                     detecta = ~t + detecta ;
                 end

                 L = L + 1 ;
                 Pcal0=eval(Poteliminado);
                 Jacobiano0=eval(Jacobiano);
                 DeltaP0=Pesp-Pcal0;
                 DeltaP0 = vpa(DeltaP0,6);
                 Delta0=(inv(Jacobiano0))*DeltaP0;
                 NRM = length(pot)+length(poti);
                 Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
                 Variable0= eval(Salidaeliminado)+(Delta0);
                 Variable0 = vpa(Variable0,6);
                 nrdp(L+1,1:length(Poteliminado)) =  DeltaP0;
                 nrds(L+1,1:length(Salidaeliminado)) = Variable0;

                 for k=1:1:length(EliminacionO)
                 eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
                 end


                 for k=1:1:length(EliminacionV)
                 eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
                 end


            end
            clc
            toc
            disp('Numero de iteraciones:')
            disp(L)
            disp('Delta de potencia por cada iteración:')
            disp([Entradaelminado.';vpa(nrdp)])
            disp('Variables Objetivo por cada iteración:')
            disp([Salidaeliminado.';vpa(nrds)])
 
    case 'NRDR'
            % Newthon Raphson desacoplado rapido
            tic
            nrdrp = zeros(1,length(Poteliminado));
            nrdrs = zeros(1,length(Salidaeliminado));
            Nm = Nm*V2;
            Lm = Lm*V2;
            Jm = zeros(size(Jm));
            Nm = zeros(size(Nm));
            Jacobianom = [Hm,Nm;Jm,Lm];
            Jacobiano = Jacobianom; 
            ReemplaYa = zeros(N-1,N-1); 
            ReemplaYv = zeros(M,M);
            ti = length(EliminacionO);
            to = length(EliminacionV);
            % Reemplazo de impedancias
            Angulos
            for z=1:1:ti
                
                for k=1:1:ti
                    
                    ReemplaYa(z,k) = -imag(Y(EliminacionO(z),EliminacionO(k)));
                end    
            
            end
            
            Voltajes
            for z=1:1:to
                
                for k=1:1:to
                    
                    ReemplaYv(z,k) = -imag(Y(EliminacionV(z),EliminacionV(k)));
                end    
            
            end
            
            Hm = ReemplaYa;
            Lm = ReemplaYv;
            Jacobianom = [Hm,Nm;Jm,Lm];
            Jacobiano = Jacobianom;
            
            Pcal0=eval(Poteliminado);
            Jacobiano0=(Jacobiano);
            DeltaP0=Pesp-Pcal0;
            DeltaP0=DeltaP0/V2;
            Delta0=(inv(Jacobiano0))*DeltaP0;
            NRM = length(pot)+length(poti);
            Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
            Variable0= eval(Salidaeliminado)+ Delta0;
            
            nrdrp(1,1:length(Poteliminado)) =  DeltaP0;
            nrdrs(1,1:length(Salidaeliminado)) = Variable0;
            
            e=0.001;
            a=600; 
            L=0; % Marca la iteracion 0
            
            for k=1:1:length(EliminacionO)
            eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
            end
            
            
            for k=1:1:length(EliminacionV)
            eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
            end
            
            detecta = 0;
            
            while Ni ~= detecta
                 detecta = 0;
                 for f = 1:1:Ni
                     t = e <= abs(DeltaP0(f,1));
                     t = eval(t);
                     detecta = ~t + detecta ;
                 end
            
                 L = L + 1 ;
                 Pcal0=eval(Poteliminado);
                 Jacobiano0=(Jacobiano)
                 DeltaP0=Pesp-Pcal0;
                 DeltaP0=DeltaP0/V2;
                 DeltaP0 = vpa(DeltaP0,6); 
                 Delta0=(inv(Jacobiano0))*DeltaP0;
                 NRM = length(pot)+length(poti);
                 Delta0(length(pot)+ 1:NRM) = (Delta0(length(pot)+1:NRM))*V2;
                 Variable0= eval(Salidaeliminado)+(Delta0);
                 Variable0 = vpa(Variable0,6);
                 nrdrp(L+1,1:length(Poteliminado)) =  DeltaP0;
                 nrdrs(L+1,1:length(Salidaeliminado)) = Variable0;
               
                 for k=1:1:length(EliminacionO)
                 eval(sprintf('O%d = Variable0(k,1) ', EliminacionO(k)));
                 end
            
            
                 for k=1:1:length(EliminacionV)
                 eval(sprintf('V%d = Variable0(k+b,1) ', EliminacionV(k)));
                 end
            
                 
            end
            clc
            toc
            disp('Numero de iteraciones:')
            disp(L)
            disp('Delta de potencia por cada iteración:')
            disp([Entradaelminado.';vpa(nrdrp)])
            disp('Variables Objetivo por cada iteración:')
            disp([Salidaeliminado.';vpa(nrdrs)])
            
    otherwise
             disp('ERROR Metodo no reconocido')
end

