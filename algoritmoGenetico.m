function [a,pie_chat_IM,t10]=algoritmoGenetico(~,areaLibre)
infanteria=true;

if infanteria
    f = waitbar(0,'Está cargando el programa, por favor espere...','Name','Estado del programa');
    
    [areaL,inverter,panel,turbina,battery,lco,clima,potencia_requerida,diesel]=cargaExcel();
    battery.eficiencia=battery.eficiencia./100;
    turbina.eficiencia=turbina.eficiencia./100;
    inverter.eficiencia=inverter.eficiencia./100;
%    diesel.eficiencia=diesel.eficiencia./100;
    rng(1,'philox');

    pot_panel=@(irradiancia,area,eficiencia_panel)(irradiancia.*eficiencia_panel.*area);
    maxPanel=@(potencia_requerida,irradiancia,area,eficiencia_panel)(ceil(potencia_requerida./(pot_panel(irradiancia,area,eficiencia_panel).*10^-3)));
    %pot_turbina=@(d_aire,area_barrido,eficiencia_turbina,vel_viento)(1/2.*d_aire.*area_barrido.*eficiencia_turbina.*vel_viento.^3);     
    pot_turbina=specialTurbine(clima,turbina);
    maxTurbina=@(potenciaNecesaria,d_aire,area_barrido,eficiencia_turbina,vel_viento)...
         (ceil(potenciaNecesaria./(pot_turbina.*10^-3)));
    rangos=[1129,1152;2113,2137;5953,5976;7369,7391]; 
    helpRandi=@(cantidadEquipo,minEquip,maxEquip)(randi([minEquip,maxEquip],cantidadEquipo,1));
    pError=@(aprox,exacto)(abs(aprox-exacto).*100./exacto);
    v_h=@(h,h_ref,v_href,alpha)((h/h_ref).^alpha.*v_href);
    %%
    %Especificaciones Algoritmo Genetico
    max_gen=500;  
    number_equip=100;
    pos_min=8.5e-4;       
    cutting=round(number_equip*0.4/2);
    prob_mutation=1/100;  
    k_friends=round(number_equip.*0.3);
    ver_24=false;
    trp=false;
    fast_mode=true;
    %%
    
    %potencia_requerida=curva_de_carga;
    clima.altura=turbina.alturaReferencia;
    h_ref=turbina.alturaReferencia;
    h=turbina.alturaUsada;
    battery.valueSelected=0;
    if h ==0
        h=10;
    end
    %areaL=areaLibre;
    alpha=turbina.alpha;
    clima.velViento=v_h(h,h_ref,clima.velViento,alpha);
    lco.total=zeros(number_equip,1);
    memory_SOCi=zeros(number_equip,1);
    memory_SOCL=zeros(number_equip,length(potencia_requerida));
    cargaPromedioBaterias=zeros(length(potencia_requerida),1);
    mP=maxPanel(potencia_requerida,clima.irradiancia,panel.area,panel.eficiencia);
    max_panel_area=ceil(areaL/panel.area);
    mP(isinf(mP))=0;
    mP(isnan(mP))=0;
    mP(mP>max_panel_area)=max_panel_area; 
    max_turbina=ceil(areaL/(turbina.areaOcupada));
    mT=maxTurbina(potencia_requerida,clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento);
    mT(isinf(mT))=0;
    mT(isnan(mT))=0; 
    mT(mT>max_turbina)=max_turbina; 
    potenciaSuplida=zeros(length(potencia_requerida),1);
    energy_accumulator=zeros(length(potencia_requerida),3);
    timer_keep=zeros(length(potencia_requerida),1);
    config=energy_accumulator;
    structure_memory=struct();
    distribucionHoras=round(linspace(1,length(potencia_requerida),ceil(length(potencia_requerida)*0.001)+1));
    valorDiesel=diesel.potencia ;
    for hora=1:length(potencia_requerida)
        tic;
        if hora ==1

            waitbar(hora/length(potencia_requerida),f,'El programa está optimizando');
        else
            timeOut=mean(timer_keep(timer_keep>0))+std(timer_keep(timer_keep>0))/2;
            timeUsed=(timeOut)*(length(potencia_requerida)-hora);
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
            %% Iteraciónes
            if generacion==1
                panel.cantidad=helpRandi(number_equip,0,mP(hora));
                panel.cantidad(end)=0;
                turbina.cantidad=helpRandi(number_equip,0,mT(hora));
                turbina.cantidad(end)=0; 
            end
            maxiCantidadT=ceil((areaL-panel.cantidad*panel.area)/(turbina.areaOcupada));
            maxiCantidadT(maxiCantidadT<0)=0;
            for i=1:length(turbina.cantidad)
                turbina.cantidad(i)=helpRandi(1,0,maxiCantidadT(i));
            end
            turbina.cantidad(end)=0;
            battery.SOCi=memory_SOCi;
            battery.SOCL=memory_SOCL;
            
            [panel,turbina,battery,diesel,lco,potenciaUsada]=planta_new(clima,panel,turbina,inverter,battery,lco,potencia_requerida,hora,diesel);
            
            
            
            memoria_lcoe(generacion)=mean(lco.total);
            memoria_equipos(generacion,:)=mean([panel.cantidad,turbina.cantidad]);
            [lco_minimo,index_mLcoe]=min(lco.total);
            best_lcoe(generacion)=lco_minimo;
            best_equipos(generacion,:)=[panel.cantidad(index_mLcoe),turbina.cantidad(index_mLcoe)];
            probability=1./(1+exp(-lco.total));
            if generacion>1
                if pError(memoria_lcoe(generacion-1),memoria_lcoe(generacion))<pos_min || generacion>max_gen
                    if ver_24==true
                        hora_ver=hora;
                    else
                        hora_ver=distribucionHoras;
                    end
                    potenciaSuplida(hora)=mean(potenciaUsada.energiaGenerada);
                    memory_SOCi=repmat(mean(battery.SOCi),length(battery.SOCi),1);
                    cargaPromedioBaterias(hora)=mean(memory_SOCi);
                    memory_SOCL(:,hora)=repmat(mean(battery.SOCL(:,hora)),length(battery.SOCL(:,hora)),1);
                    config(hora,:)=mean([panel.cantidad,turbina.cantidad,lco.total]);
                    energy_accumulator(hora,:)=[mean([abs(battery.SOCL(:,hora)),diesel.generar]),generacion];
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
            poblation_mutation1=mutacion(poblation1,prob_mutation);
            poblation_mutation2=mutacion(poblation2,prob_mutation);

            panel.cantidad=bi2de(poblation_mutation1);
            turbina.cantidad=bi2de(poblation_mutation2);

            panel.cantidad(panel.cantidad>mP(hora))=mP(hora);
            %turbina.cantidad(panel.cantidad>mP(hora))=0;
            
            turbina.cantidad(turbina.cantidad>mT(hora))=mT(hora);
            %panel.cantidad(turbina.cantidad>mT(hora))=0;
            
            generacion=generacion+1;
        end

        timer_keep(hora)=toc;
    end
    %%
    figure ('Name','Cambio carga baterias')
    plot(cargaPromedioBaterias,'g-o');
    xlabel("Horas");
    ylabel("EnergiaBaterias");
    title('Carga de baterias')
    grid()
    %%
    
    %
    topPot=max(potencia_requerida);
    k1=config(:,1)>0; %Las configuraciones que tienen más de 0 paneles
    %En todos usa potencia nominal [kW]
    energySun=panel.potencia.*config(:,1); %<- Potencia generada por panel
    eneryWind= turbina.potencia.*config(:,2);%<- Potencia generada por turbina
    energyBattery=ceil(energy_accumulator(:,1)./battery.SOCMax).*battery.amperioHora.*battery.voltaje.*10^-3; %<- Potencia generada por bateria
    eneryDiesel=energy_accumulator(:,2);%<- Potencia generada por diesel
    mk=(energySun+eneryWind+energyBattery+eneryDiesel)./topPot; %Compara todas las energias con la potencia máxima
    
    [~,idS]=sort(logical(k1.*(mk<=1)),'descend');
    idS=idS(1);
    %a=[paneles.cantidad(id),turbina.cantidad(id),energy_accumulator(id,1:2),horas_activo,median(energy_accumulator(:,3))];
    switch trp
        case false

            horas_activo=sum(energy_accumulator(:,2)>0);
            %a=[ceil(paneles_cantidad),ceil(mean(config(:,2))),mean(config(:,3)),ceil(sum(energy_accumulator(:,1),1)./battery.SOCMax),ceil(max(energy_accumulator(:,2))./valorDiesel),horas_activo,median(energy_accumulator(:,3))];
            a=[config(idS,:),sum(energy_accumulator(:,1:2)),horas_activo,median(energy_accumulator(:,3))];
            battery.valueSelected=ceil(sum(energy_accumulator(:,1),1)./battery.SOCMax);
            disp(ceil(sum(energy_accumulator(:,2),1)));
        case true
            a=[config(hora_ver,:),...
                ceil(sum(energy_accumulator(:,1),1)./battery.SOCMax),sum(energy_accumulator(:,2),1),horas_activo,median(energy_accumulator(:,3))];
    end

    a=array2table(a,'variableNames',{'Modulo','Turbina','LCOE','Numero bateria','Numero motores diesel','Horas diesel activo','IteracionMediana'});
    a.Modulo=ceil(a.Modulo);
    a.Turbina=ceil(a.Turbina);
    a.("Numero motores diesel")=ceil(a.("Numero motores diesel")./valorDiesel.*10^-3);
    a.("Numero bateria")=ceil(a.("Numero bateria") ./(battery.SOCMax.*10^3));
%     setM=normalize([a.Modulo,a.Turbina,a.LCOE]);
%     setA=normalize(config);
%     distancesset=pdist2(setM,setA);
%     [~,horaWin]=min(distancesset);
%     a.Modulo=ceil(config(horaWin,1));
%     a.Turbina=ceil(config(horaWin,2));
%     a.LCOE=config(horaWin,3);
    horaWin=idS;
    printImages(horaWin,structure_memory(horaWin).memoria_equipos,structure_memory(horaWin).config,...
        structure_memory(horaWin).best_equipos,structure_memory(horaWin).best_lcoe,structure_memory(horaWin).memoria_lcoe);
    disp(a);
    panel.cantidad=a.Modulo;
    turbina.cantidad=a.Turbina;
    memory_SOCi=zeros(1,1);
    memory_SOCL=zeros(1,length(potencia_requerida));
    battery.SOCi=memory_SOCi;
    battery.SOCL=memory_SOCL;
    keep_ans=zeros(length(potencia_requerida),6);
    keep_ans(:,2)=potencia_requerida;
    memoria_batterias_u=[];
    battery.maxEnergySelected=battery.valueSelected*battery.SOCMax;
    maxDiesel= sum(energy_accumulator(:,2),1)*1.5;
    id=1;
    battery.carga=0;
    battery.descarga=0;
    for hora=1:length(potencia_requerida)
        
        battery.SOCi=memory_SOCi;
        [panel,turbina,battery,diesel,lco,potenciaUsada]=planta_new(clima,panel,turbina,inverter,battery,lco,potencia_requerida,hora,diesel);
        memory_SOCi=battery.SOCi;
%         maxDiesel=maxDiesel-potenciaUsada.diesel;
%         if maxDiesel<=0
%             potenciaUsada.diesel=0;
%             id=id+1;
%         end
        %potenciaRenovable=potenciaUsada.panel+potenciaUsada.turbina+battery.SOCi;
        potenciaGenerada=potenciaUsada.energiaGenerada;
        keep_ans(hora,2)=keep_ans(hora,2)+battery.carga; %
        assert(round(potenciaGenerada,3)==round(keep_ans(hora,2),3));
        keep_ans(hora,[1,3:end])=[potenciaGenerada,...
                                  potenciaUsada.panel,...
                                  potenciaUsada.turbina,...
                                  battery.descarga,...
                                  potenciaUsada.diesel];
        if battery.SOCi>battery.maxEnergySelected
            battery.SOCi=battery.maxEnergySelected;
        end
    end
     variablesName={'EnergiaModulo','EnergiaTurbina','EnergiaBaterias','EnergiaMotor','EnergiaRequerida','EnergiaGenerada'};
     keep_ansTable=array2table([keep_ans(:,3:end),keep_ans(:,1),keep_ans(:,end)],'VariableNames',variablesName);
     keep_ans=keep_ans(:,[3:end,2,1]);
    %%
    disp(sum(energy_accumulator(:,1)))
    pie_chat_IM=figure ('Name','Diagrama de distribución energetica');
    if true
        pPused=a.Modulo*panel.potencia ; 
        pTUsed=a.Turbina*turbina.potencia;
        pAcumUsed=max(keep_ansTable.('EnergiaMotor'));
        disp("La energía del motor es:");
        disp(pAcumUsed);
    else
        pPused=sum(keep_ansTable.('EnergiaModulo'));
        pTUsed=sum(keep_ansTable.('EnergiaTurbina'));
        pAcumUsed=max(keep_ansTable.('EnergiaMotor'));
        bAcumUsed=sum(energy_accumulator(:,1),1)*10^3;
    end
    if true
        vectorPotencias=[pPused;pTUsed;pAcumUsed];
        labelsPorcentajes={'Modulos PV ','Turbinas eolicas ','Generador(es) Diesel '}';
    else
        vectorPotencias=[pPused;pTUsed;pAcumUsed;bAcumUsed];
        labelsPorcentajes={'Modulos PV ','Turbinas eolicas ','Generador(es) Diesel ','Baterias'}';
    end

    porcentajes=vectorPotencias/sum(round(vectorPotencias,2));
    %porcentajes=porcentajes(1)+(1-sum(porcentajes));
    %%
    pie_porcentajes=pie(porcentajes);
    colormap ([1,1,0;
               0,0,1;
               0,1,0;])
    pText = findobj(pie_porcentajes,'Type','text');
    percentValues = get(pText,'String'); 

    combinedtxt = strcat(labelsPorcentajes,percentValues); 
    for i=1:length(labelsPorcentajes)
     pText(i).String=combinedtxt(i);
    end
    %%
    try
        assert(length(keep_ans)/24>1);
        packData=zeros(length(keep_ans)/24,1);
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
    tabla_energias=array2table(vPotenciasNew,'VariableNames',variablesName,'RowName',diaName);
    disp(tabla_energias)
    total=true;   
    figure('Name','Tabla de validación de la herramienta 1');
    title('Validación de la herramienta');
    if ~total
        variablesNameUse=variablesName([end-1,end]);
    else
        variablesNameUse=variablesName;
    end
    j=1;
    
    formas=["-o","-s"];
    for name = variablesNameUse(1:end-2)
        p=randi(length(formas),1);
        color='ycbg';
        plot(tabla_energias.(string(name)),formas(p),'Color',color(j),'MarkerFaceColor',color(j),...
            'MarkerEdgeColor','k','DisplayName',string(name));
        j=j+1;
        hold on
    end
    clear j
    grid minor ;
    legend();
    xlabel('Día');
    ylabel('Energía kWh');
    figure('Name','Tabla de validación de la herramienta 1');
    title('Validación de la herramienta');
    j=1;
        formas=["-o","-s"];
    for name = variablesNameUse(end-1:end)
        p=randi(length(formas),1);
        color='rb';
        plot(tabla_energias.(string(name)),formas(p),'Color',color(j),'MarkerFaceColor',color(j),...
            'MarkerEdgeColor','k','DisplayName',string(name));
        j=j+1;
        hold on
    end
    clear j
    grid minor ;
    legend();
    xlabel('Día');
    ylabel('Energía kWh');
 
    [x_rango,~]=size(rangos);
    diasMuestra=zeros(x_rango,2);
    diaName=string(zeros(x_rango,1));
    try
        for i=1:x_rango 
            diaName(i)="Día "+string(i);
            rango_1=rangos(i,1);
            rango_2=rangos(i,2);
            diasMuestra(i,:)=[sum(keep_ans(rango_1:rango_2,end-1),1),vPotenciasNew(:,end)];
        end
        tabla_muestra=array2table(diasMuestra,'VariableNames',variablesName(end-1:end),'RowName',diaName);
        disp(tabla_muestra);
    catch
        warning("Se requieren más datos");
    end


    try
    %%
    tiempo=1;
    panel.energiaTiempo=sum(pot_panel(clima.irradiancia,panel.area,panel.eficiencia)).*a.Modulo.*tiempo.*10^-3+sum(energy_accumulator(:,1)); %kWh/mes
    turbina.energiaTiempo=sum(pot_turbina).*a.Turbina.*tiempo.*10^-3;
    %turbina.energiaTiempo=sum(pot_turbina(clima.densidadAire,turbina.areaBarrido,turbina.eficiencia,clima.velViento)).*a.Turbina.*tiempo.*10^-3;
    dieselU.energiaTiempo=pAcumUsed(:,1).*10^-3.*tiempo;
    t10=table(panel.energiaTiempo,turbina.energiaTiempo,dieselU.energiaTiempo,'VariableNames',{'PVkWhmonth','WindkWhmonth','dieselkWhmonth'});
    close (f);

    catch ME
        close(f);
        messageError=sprintf("Ha ocurrido un error\n %s",ME.message);
        errordlg(messageError,'Error');
    end
end
%Cruce y mutación
end
function [generacionNew]=reproduccion(genetic,results,k_friends,k_child)
[row,col,deep]=size(genetic);
incubator=zeros(k_child*2,col,deep);
generacionNew=zeros(row,col,deep);
if k_friends>row
    k_friends=row-1;
end
ventaja=round(k_child*0.01);
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
        if length(puntuation)==1
            puntuation=results(site);
        end
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
function [pT]=specialTurbine(clima,turbina)
    pT=clima.velViento*0;
    c1=turbina.velocidadArranque<=clima.velViento;
    c2=turbina.velocidadNominal>clima.velViento ;
    c3=clima.velViento<=turbina.velocidadMaxima ;
    v=clima.velViento(logical(c1.*c2));
    pT(logical(c1.*c2))=v.^3 .*(turbina.potencia./(turbina.velocidadNominal.^3-turbina.velocidadArranque.^3))-turbina.potencia.*(turbina.velocidadArranque.^3./(turbina.velocidadNominal.^3-turbina.velocidadArranque.^3));
    pT(logical(~c2.*c3))=turbina.potencia;
    pT=pT.*turbina.eficiencia;
    return;
end