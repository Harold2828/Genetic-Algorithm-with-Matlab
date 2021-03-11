function [a,pie_chat_IM,t10]=algoritmoGenetico(curva_de_carga,altura,areaLibre)
f = waitbar(0,'Está cargando el programa, por favor espere...');
try
rng(1);
%Ecuaciones
% maxPanel=@(eficiencia,areaLibre,irradiancia,potenciaPanel)...
%     (eficiencia.*areaLibre.*irradiancia.*10^-3./potenciaPanel);
%Para la potencia del panel
pot_panel=@(irradiancia,area,eficiencia_panel)(irradiancia.*eficiencia_panel.*area);
%Maxima cantidad de paneles
maxPanel=@(potencia_requerida,irradiancia,area,eficiencia_panel)(ceil(potencia_requerida./(pot_panel(irradiancia,area,eficiencia_panel).*10^-3)));
%Maxima turbina
%Para la potencia de la turbina
 pot_turbina=@(d_aire,area_barrido,eficiencia_turbina,vel_viento)(1/2.*d_aire.*area_barrido.*eficiencia_turbina.*vel_viento.^3);
%Para la maxima cantidad de turbinas     
 maxTurbina=@(potenciaNecesaria,d_aire,area_barrido,eficiencia_turbina,vel_viento)...
     (ceil(potenciaNecesaria./(pot_turbina(d_aire,area_barrido,eficiencia_turbina,vel_viento).*10^-3)));
rangos=[1129,1152;2113,2137;5953,5976;7369,7391];
helpRandi=@(cantidadEquipo,minEquip,maxEquip)(randi([minEquip,maxEquip],cantidadEquipo,1));
%Porcentaje de error
pError=@(aprox,exacto)(abs(aprox-exacto).*100./exacto);
%Velocidad del viento
v_h=@(h,h_ref,v_href,alpha)((h/h_ref).^alpha.*v_href);
%%
%Especificaciones Algoritmo Genetico
max_gen=200; %Por qué 200 ?
pos_min=1e-5;
number_equip=700;
cutting=round(number_equip*0.35/2);
prob_mutation=0.01;
k_friends=round(number_equip.*0.35);
ver_24=false;
trp=false;
fast_mode=true;
%%
[~,inverter,panel,turbina,battery,lco,clima,~]=cargaExcel();
potencia_requerida=curva_de_carga;
clima.altura=turbina.alturaReferencia;
h_ref=turbina.alturaReferencia;
h=turbina.alturaUsada;
if h ==0
    h=10;
end
areaL=areaLibre;

alpha=turbina.alpha;
clima.velViento=v_h(h,h_ref,clima.velViento,alpha);
lco.total=zeros(number_equip,1);
memory_SOCi=zeros(number_equip,1);
memory_SOCL=zeros(number_equip,length(potencia_requerida));

pot_tubina_1=pot_turbina(clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento).*10^-3;
pot_panel_1=pot_panel(clima.irradiancia,panel.area,panel.eficiencia).*10^-3;
potencia_sePuedeGenerar=pot_tubina_1+pot_panel_1;
potencia_requerida(potencia_requerida<=0)=potencia_sePuedeGenerar(potencia_requerida<=0);
%Cantidades maximas
mP=maxPanel(potencia_requerida,clima.irradiancia,panel.area,panel.eficiencia);
max_panel_area=ceil(areaL/panel.area);
mP(isinf(mP))=0;
mP(isnan(mP))=0;
mP(mP>max_panel_area)=max_panel_area; 
%Maxima turbina
max_turbina=ceil(areaL/(153 * 27*0.01));
mT=maxTurbina(potencia_requerida,clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento);
mT(isinf(mT))=0;
mT(isnan(mT))=0; 
mT(mT>max_turbina)=max_turbina; 
potenciaSuplida=zeros(length(potencia_requerida),1);
energy_accumulator=zeros(length(potencia_requerida),3);
config=energy_accumulator;
structure_memory=struct();
distribucionHoras=round(linspace(1,length(potencia_requerida),ceil(length(potencia_requerida)*0.001)+1));
for hora=1:length(potencia_requerida)
    if hora ==1
        tic;
        waitbar(hora/length(potencia_requerida),f,'El programa está optimizando');
    else
        timeUsed=(timeOut+2)*length(potencia_requerida);
        message2='segundos';
        if timeUsed==1
            message2='segundo';
        end
        if timeUsed>60
            timeUsed=timeUsed/60;
            message2='minutos';
            if timeUsed>60
                timeUsed=timeUsed/60;
                message2='horas';
            end
        end
        messageWaitbar=sprintf('El programa está optimizando\nEl tiempo de espera estimado es %d %s',char(timeUsed),message2);
        waitbar(hora/length(potencia_requerida),f,messageWaitbar);
    end
    
    
    generacion=1;
    memoria_lcoe=zeros(1,1);
    memoria_equipos=zeros(1,2);
    best_lcoe=zeros(1,1);
    best_equipos=zeros(1,2);
    if logical(sum(hora==distribucionHoras))
    end
    while true
        %%
        if generacion==1
            
            %Combinaciones de equipos
            panel.cantidad=helpRandi(number_equip,0,mP(hora));
            panel.cantidad(end)=0;
            turbina.cantidad=helpRandi(number_equip,0,mT(hora));
            turbina.cantidad(end)=0; 
        end
        battery.SOCi=memory_SOCi;
        battery.SOCL=memory_SOCL;
        [panel,turbina,battery,diesel,lco,potenciaUsada]=planta_new(clima,panel,turbina,inverter,battery,lco,potencia_requerida,hora);
        memoria_lcoe(generacion)=mean(lco.total);
        memoria_equipos(generacion,:)=mean([panel.cantidad,turbina.cantidad]);
        [lco_minimo,index_mLcoe]=min(lco.total);
        best_lcoe(generacion)=lco_minimo;
        best_equipos(generacion,:)=[panel.cantidad(index_mLcoe),turbina.cantidad(index_mLcoe)];
        probability=log(lco.total);
        if generacion>1
            if pError(memoria_lcoe(generacion-1),memoria_lcoe(generacion))<pos_min || generacion>max_gen
                if ver_24==true
                    hora_ver=hora;
                else
                    hora_ver=distribucionHoras;
                end
                potenciaSuplida(hora)=mean(potenciaUsada.energiaGenerada);
                memory_SOCi=ones(number_equip,1).*mean(battery.SOCi(~isoutlier(battery.SOCi)));
                mem_SOCLH=ones(number_equip,1).*mean(battery.SOCL(~isoutlier(battery.SOCL(:,hora)),hora));
                memory_SOCL(:,hora)=mem_SOCLH;
                
                config(hora,:)=mean([panel.cantidad,turbina.cantidad,lco.total]);
                energy_accumulator(hora,:)=[mean([abs(battery.SOCL(:,hora)),diesel]),generacion];
                %Incluyendo diesel 100%
                structure_memory(hora).memoria_equipos=memoria_equipos;
                structure_memory(hora).config=config;
                structure_memory(hora).best_equipos=best_equipos;
                structure_memory(hora).best_lcoe=best_lcoe;
                structure_memory(hora).memoria_lcoe=memoria_lcoe;
                if ~fast_mode
                    printImages(hora,memoria_equipos,config,best_equipos,best_lcoe,memoria_lcoe)
                end
                break;
            end
        end
        if sum(panel.cantidad==0)==length(panel.cantidad)
            part1=panel.cantidad;
        else
            part1=de2bi(panel.cantidad);
        end
        
        if sum(turbina.cantidad==0)==length(turbina.cantidad)
            part2=turbina.cantidad;
        else
            part2=de2bi(turbina.cantidad);
        end
           
        %Cruce
        k_child=cutting;
        poblation1=reproduccion(part1,probability,k_friends,k_child);
        poblation2=reproduccion(part2,probability,k_friends,k_child);
        %Mutacion
        poblation_mutation1=mutacion(poblation1,prob_mutation);
        poblation_mutation2=mutacion(poblation2,prob_mutation);

        panel.cantidad=bi2de(poblation_mutation1);
        turbina.cantidad=bi2de(poblation_mutation2);
        
        panel.cantidad(panel.cantidad>mP(hora))=mP(hora);
        turbina.cantidad(turbina.cantidad>mT(hora))=mT(hora);
        generacion=generacion+1;
    end
    if hora==1
        timeOut=toc;
    end
end
%Aquí finaliza el procedimiento
%Comienza el resto, tabalas.. imagenes...

switch trp
    case false

        horas_activo=sum(energy_accumulator(:,2)>0);
        if isnan(mean(config(config(:,1)>0,1)))
            paneles_cantidad=0;
        else
            paneles_cantidad=mean(config(config(:,1)>0,1));
        end
        a=[ceil(paneles_cantidad),ceil(mean(config(:,2))),mean(config(:,3)),ceil(sum(energy_accumulator(:,1),1)./battery.SOCMax),sum(energy_accumulator(:,2),1),horas_activo,median(energy_accumulator(:,3))];
    
    case true
        a=[config(hora_ver,:),...
            ceil(sum(energy_accumulator(:,1),1)./battery.SOCMax),sum(energy_accumulator(:,2),1),horas_activo,median(energy_accumulator(:,3))];
end
a=array2table(a,'variableNames',{'Panel','Turbina','LCOE','Numero bateria','Energia diesel','Horas diesel activo','IteracionMediana'});
for hora=1:length(structure_memory)
    operatorMin=pError(mean([structure_memory(hora).memoria_equipos(end,:),...
        structure_memory(hora).memoria_lcoe(end)]),mean([a.Panel,a.Turbina,a.LCOE]));
    if hora==1
        horaWin=hora;
        minValue=operatorMin;
    elseif minValue>operatorMin
        horaWin=hora;
        minValue=operatorMin;
        
    end
end
%printImages(hora,memoria_equipos,config,best_equipos,best_lcoe,memoria_lcoe)
printImages(horaWin,structure_memory(horaWin).memoria_equipos,structure_memory(horaWin).config,...
    structure_memory(horaWin).best_equipos,structure_memory(horaWin).best_lcoe,structure_memory(horaWin).memoria_lcoe);
disp(a);
panel.cantidad=a.Panel;
turbina.cantidad=a.Turbina;
memory_SOCi=zeros(1,1);
memory_SOCL=zeros(1,length(potencia_requerida));
battery.SOCi=memory_SOCi;
battery.SOCL=memory_SOCL;
keep_ans=zeros(length(potencia_requerida),6);
keep_ans(:,2)=potencia_requerida;

for hora=1:length(potencia_requerida)
    [panel,turbina,battery,diesel,lco,potenciaUsada]=planta_new(clima,panel,turbina,inverter,battery,lco,potencia_requerida,hora);
    keep_ans(hora,[1,3:end])=[potenciaUsada.energiaGenerada,...
                              potenciaUsada.panel,...
                              potenciaUsada.turbina,...
                              battery.SOCL(hora),...
                              potenciaUsada.diesel];
end
 keep_ansTable=array2table(keep_ans(:,1:2),'VariableNames',{'EnergiaGenerada','EnergiaRequerida'});
 disp(keep_ansTable); %La tabla de generación de energia con la config. seleccionada
 keep_ans=keep_ans(:,[3:end,2,1]);
%%
if ver_24
    [~,idx_sort]=sort(pdist2([a.Panel,a.Turbina,a.LCOE],config));
    keep_figures=[idx_sort(1)*2-1,idx_sort(1)*2];
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, keep_figures)); %
end

disp(sum(energy_accumulator(:,1)))
pie_chat_IM=figure ('Name','Pie chart de energia');
pPused=sum(pot_panel(clima.irradiancia,panel.area,panel.eficiencia).*config(:,1));
pTUsed=sum(pot_turbina(clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento).*config(:,2));
pAcumUsed=sum(energy_accumulator(:,2),1)*10^3;
vectorPotencias=[pPused;pTUsed;pAcumUsed'];
porcentajes=vectorPotencias/sum(vectorPotencias);
%%
pie_porcentajes=pie(porcentajes);
pText = findobj(pie_porcentajes,'Type','text');
percentValues = get(pText,'String'); 
labelsPorcentajes={'Paneles ','Turbinas eolicas ','Generador Diesel '}';
combinedtxt = strcat(labelsPorcentajes,percentValues); 
for i=1:length(labelsPorcentajes)
 pText(i).String=combinedtxt(i);
end
%%

%vPotencias=[pot_panel(clima.irradiancia,panel.area,panel.eficiencia).*config(:,1).*10^-3,...
%    pot_turbina(clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento).*config(:,2).*10^-3,...
%    energy_accumulator(:,1:2),potencia_requerida];

try
    i=1;
    while true
        pack1=24*(i-1);
        if pack1==0
            pack1=1;
            vPotenciasNew=[];
        end
        pack2=24*i;
        if i<=(length(keep_ans)/24)
            vPotenciasNew=[vPotenciasNew;sum(keep_ans(pack1:pack2,:),1)];

            i=i+1;
        else
            i=i-1;
            diaName=repmat("Dia ",i,1)+string(1:i)';
            break;
        end
        
    end
catch
    disp("Error en las imagenes, se mostraran todas las horas");
    vPotenciasNew=keep_ans;
    diaName=repmat("Hora ",length(keep_ans),1)+string((1:length(keep_ans)))';
end
vPotenciasNew(:,2)=repmat(mean(vPotenciasNew(:,2)),length(vPotenciasNew(:,2)),1);
variablesName={'EnergiaPanel','EnergiaTurbina','EnergiaBaterias','EnergiaMotor','EnergiaRequerida','EnergiaGenerada'};
tabla_energias=array2table(vPotenciasNew,'VariableNames',variablesName,'RowName',diaName);
disp(tabla_energias)

[x_rango,~]=size(rangos);
diasMuestra=zeros(x_rango,2);
diaName=string(zeros(x_rango,1));
try
    for i=1:x_rango %Tiene que estar presente los 8000 datos
        diaName(i)="Dia "+string(i);
        rango_1=rangos(i,1);
        rango_2=rangos(i,2);
        diasMuestra(i,:)=sum(keep_ans(rango_1:rango_2,end-1:end),1);
    end
    tabla_muestra=array2table(diasMuestra,'VariableNames',variablesName(end-1:end),'RowName',diaName);
    disp(tabla_muestra);
catch
    warning("Se requieren más datos");
end
    %La grafica de las líneas
    figure('Name','Tabla de Comparación de las energías');
    title('Comparación de energía');
    hold on
    plot(tabla_energias.EnergiaRequerida,'b-o','DisplayName','EnergiaRequerida');
    plot(tabla_energias.EnergiaGenerada,'k--','DisplayName','EnergiaGenerada');
    grid;
    legend();
    xlabel('Hora');
    ylabel('Energía kW');



%%
%para las tablas
tiempo=1;
panel.energiaTiempo=sum(pot_panel(clima.irradiancia,panel.area,panel.eficiencia)).*a.Panel.*tiempo.*10^-3+sum(energy_accumulator(:,1)); %kWh/mes
turbina.energiaTiempo=sum(pot_turbina(clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento)).*a.Turbina.*tiempo.*10^-3;
dieselU.energiaTiempo=pAcumUsed(:,1).*10^-3.*tiempo;

t10=table(panel.energiaTiempo,turbina.energiaTiempo,dieselU.energiaTiempo,'VariableNames',{'PVkWhmonth','WindkWhmonth','dieselkWhmonth'});
%Aquií termina
close (f);
catch ME
    close(f);

    messageError=sprintf("Ha ocurrido un error\n %s",ME.message);
    f=errordlg(messageError,'Error');
    pause(5)
    close(f)
end
end
%Para el cruce
function [generacionNew]=reproduccion(genetic,results,k_friends,k_child)
[row,col,deep]=size(genetic);
incubator=zeros(k_child*2,col,deep);
generacionNew=zeros(row,col,deep);
if k_friends>row
    k_friends=row-1;
end
ventaja=round(k_child*0.21);
engagement=randi([1,row],k_friends,k_child);
[~,idx_min_LCOE]=min(results);
coord_new_lcoe=[randperm(k_friends,ventaja)',randperm(k_child,ventaja)'];
for i=1:ventaja
    coord_row=coord_new_lcoe(i,1);
    coord_col=coord_new_lcoe(i,2);
    engagement(coord_row,coord_col)=idx_min_LCOE;
end

for dim=1:deep
    st1=1;
    st2=2;
    for child=1:k_child
        site=engagement(:,child);
        puntuation=results(unique(site));
        [~,idx_puntuation]=sort(puntuation);
        parents=idx_puntuation(1:2);
        armada=genetic(parents,:,dim);
        corte_gen=randi([1,col],1,1);
        incubator(st1:st2,:,dim)=[armada(1,1:corte_gen),armada(2,corte_gen+1:end);...
            armada(2,1:corte_gen),armada(2,corte_gen+1:end)];
        st1=st2+1;
        st2=st2+2;
    end
    
    [~,idx_results_sort]=sort(results);
    points=row-k_child*2;
    if points>0
        old_gen=genetic(idx_results_sort(1:row-k_child*2),:,dim);
        incubator_send=incubator(:,:,dim);
    elseif points==0
        old_gen=[];
        incubator_send=incubator(:,:,dim);
    elseif points<0
        old_gen=[];
            incubator_send=incubator(logical(mod((1:k_child*2)',2)),:,dim);
    end
    generacionNew(:,:,dim)=[old_gen;incubator_send];
end
end
%%
%Para la mutacion
function [zerg]=mutacion(genetic,probability_mutation)
genetic=logical(genetic);
if probability_mutation<1
    probability_mutation=abs(probability_mutation*100);
elseif probability_mutation>100
    probability_mutation=100;
end
[row,col,deep]=size(genetic);
number_random=randi(100,row,1,deep);
for dim=1:deep
    point=(number_random(:,:,dim)<=probability_mutation);
    selected=genetic(point,:,dim);
    [row_selected,~]=size(selected);
    quantity_mutation=randi([0,col],row_selected,1);
    for i=1:row_selected
        point_mutation=randi(col,1,unique(quantity_mutation(i)));
        selected(i,point_mutation)=~selected(i,point_mutation);
    end
    genetic(point,:,dim)=selected;
end
zerg=genetic;
end
function printImages(hora,memoria_equipos,config,best_equipos,best_lcoe,memoria_lcoe)
    mensajeUsing=sprintf("Graficas LCOE para la hora %d",hora);
    figure ('Name',mensajeUsing)
    plot3(memoria_equipos(:,2),memoria_equipos(:,1),memoria_lcoe,'k-')
    xlabel("Turbinas")
    ylabel("Panel")
    zlabel("LCOE")
    grid
    hold on
    plot3(config(hora,2),config(hora,1),config(hora,3),'go','MarkerFaceColor','g')
    plot3(best_equipos(:,2),best_equipos(:,1),best_lcoe,'r+')
    legend('Promedio del grupo','Configuracion seleccionada','Mejor individuo')
    figure ("Name","Promedio Vs Mejor individuo") %Por ahora dejemos esto asÃ­
    plot(best_lcoe,'b-o')
    hold on
    plot(memoria_lcoe,'ro','MarkerFaceColor','r')
    grid
    xlabel("Iteracion")
    ylabel("LCOE")
    legend('Mejor LCOE','LCOE promedio')
end
