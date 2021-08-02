function [panel,turbina,battery,diesel,lco,potencia]=planta_new(clima,panel,turbina,inverter,battery,lco,potenciaRequeria,hora,diesel)
%%
p_Turbina=@(eficiencia_turbina,d_aire,area_barrido,vel_viento,numero_turbina)...
    (1/2.*d_aire.*area_barrido.*eficiencia_turbina.*(vel_viento.^3).*numero_turbina.*10^-3 );

p_Panel=@(eficiencia_panel,area,irradiancia,numero_panel)...
     (irradiancia.*eficiencia_panel.*area.*numero_panel.*10^-3);

p_Diesel=@(potencia_requerida,pot_renovable)(potencia_requerida-pot_renovable);

lcoe=@(lcoSol,lcoViento,lcoBat,lcoDiesel,energy_sol,energy_wind,energy_bat,energy_diesel)((lcoSol+lcoViento+lcoBat+lcoDiesel)./eGen);
lcoe_formula=@(lcoe_used,energy_used)(lcoe_used.*energy_used);
%LCOE
I=@(costoEquipo,numeroEquipos)(costoEquipo.*numeroEquipos);
ACi=@(inversion,tasaInteres,vidaUtil)( (tasaInteres.*inversion.*(1+tasaInteres).^vidaUtil)./((1+tasaInteres).^vidaUtil -1) );
AVs=@(valorSalvamento,tasaInteres,vidaUtil)( (valorSalvamento.*tasaInteres)./((1+tasaInteres).^vidaUtil-1) );
ACcomb=@(energiaGenerador,eficienciaGenerador,heatRate,costoDiesel)( (energiaGenerador.*heatRate.*costoDiesel)./eficienciaGenerador );
%%
panel_run=panel.cantidad;
turbina_run=turbina.cantidad;
pv_gen=p_Panel(panel.eficiencia(hora),panel.area,clima.irradiancia(hora),panel_run);
tur_gen=p_Turbina(turbina.eficiencia,clima.densidadAire(hora),clima.velViento(hora),turbina_run,turbina.areaBarrido);
renovable_gen=sum([pv_gen,tur_gen],2);
posGG=logical(renovable_gen>potenciaRequeria(hora));
posUG=~posGG;
%%
if isempty(battery.SOCi(posGG))==0
    posMM_1=ones(length(battery.SOCi),1); %
    posMM=logical(posMM_1.*posGG);
    if isempty(battery.SOCi(posMM))==0
        battery.SOCi(posMM)=bateria(inverter.eficiencia,battery.eficiencia,...
            battery.SOCi(posMM),battery.SOCMax,battery.SOCMin,battery.autoDescarga,renovable_gen(posMM),...
            potenciaRequeria(hora),"carga");
        battery.carga=battery.SOCi;
        battery.carga(~posMM)=0;
    end
end
%%
if isempty(battery.SOCi(posUG))==0
    posMm_1=battery.SOCi>battery.SOCMin;
    posMm=logical(posMm_1.*posUG); 
    if isempty(battery.SOCi(posMm))==0
        antes=battery.SOCi(posMm);
        battery.SOCi(posMm)=bateria(inverter.eficiencia,battery.eficiencia,battery.SOCi(posMm),battery.SOCMax,...
            battery.SOCMin,battery.autoDescarga,renovable_gen(posMm),potenciaRequeria(hora),"descarga");
         battery.SOCL(posMm,hora)=antes-battery.SOCi(posMm);
        renovable_gen(posMm,:)=renovable_gen(posMm,:)+battery.SOCi(posMm);
        battery.descarga=battery.SOCi;
        battery.descarga(~posMm)=0;
    end
end
%%
diesel.generar=p_Diesel(potenciaRequeria(hora),renovable_gen);
diesel.generar(diesel.generar<0)=diesel.generar(diesel.generar<0).*0;
diesel.eficiencia=diesel.generar./(diesel.consumoCalorifico.*10^3);
energia_generada=renovable_gen+diesel.generar;
potencia.panel=pv_gen;
potencia.turbina=tur_gen;
potencia.diesel=diesel.generar;
%inicioLCOE
lcoeValue=zeros(length(panel.cantidad),4);
%Panel
lcoeValue(:,1)=ACi(I(panel.costo,panel.cantidad),8/100,panel.vidaUtil)+...
    I(panel.costo,panel.cantidad).*panel.coym;
%lcoeValue(:,1)=lcoeValue(:,1)./pv_gen;
%Turbina
lcoeValue(:,2)=ACi(I(turbina.costo,turbina.cantidad),8/100,turbina.vidaUtil)+...
    I(turbina.costo,turbina.cantidad).*turbina.coym;
%lcoeValue(:,2)=lcoeValue(:,2)./tur_gen;
%Bateria
lcoeValue(:,3)=ACi(I(battery.costo,ceil(battery.SOCL(:,hora)./(battery.SOCMax.*10^3))),8/100,battery.vidaUtil)+...
    I(battery.costo,ceil(battery.SOCL(:,hora)./(battery.SOCMax.*10^3))).*battery.coym+ceil(battery.SOCL(:,hora)./battery.SOCMax.*10^3).*battery.costo;
%lcoeValue(:,3)=lcoeValue(:,3)./battery.SOCL(:,hora);
%Diesel
lcoeValue(:,4)=ACi(I(diesel.costo,ceil(diesel.generar./(diesel.potencia.*10^3))),8/100,diesel.vidaUtil)+...
    I(diesel.costo,ceil(diesel.generar./(diesel.potencia.*10^3))).*diesel.coym+...
    ACcomb(diesel.generar,diesel.eficiencia,diesel.hr,diesel.costo)-...
    AVs(diesel.valorSalvamento,I(diesel.costo,ceil(diesel.generar./(diesel.potencia.*10^3))),diesel.vidaUtil);

lcoeValue(isnan(lcoeValue))=0;
lcoeValue(isinf(lcoeValue))=0;
%lcoeValue(:,4)=lcoeValue(:,4)./diesel.generar;
%Fin
potencias=[pv_gen,tur_gen,battery.SOCi,diesel.generar];
lcoe_using=[lco.sol,lco.viento,lco.bat,lco.diesel];
[col_p,~]=size(pv_gen);
lco_temporal=zeros(col_p,4);

lco.total=sum(lcoeValue,2)./energia_generada;
lco.total(isnan(lco.total))=0;
potencia.energiaGenerada=energia_generada;
try
    round_please=1;
    rememberEnergy=sum(~(round(energia_generada,round_please)>=round(potenciaRequeria(hora),round_please)));
    assert (rememberEnergy==0);
catch
end
end
 function eBateria=bateria(ef_inverter,ef_battery,SOC,SOCMax,SOCMin,ind_auto,e_gen,e_demand,tipo)
 carga=@(SOC,ind_auto,e_gen,e_demand,ef_inverter,ef_battery)...
        (SOC.*(1-ind_auto)+ef_battery.*((e_gen-e_demand)./ef_inverter));
    descarga=@(SOC,ind_auto,e_gen,e_demand,ef_inverter)...
        (SOC.*(1-ind_auto)+((-e_gen+e_demand)./ef_inverter));
    switch string(upper(tipo))
        case "CARGA"
            eBateria=carga(SOC,ind_auto,e_gen,e_demand,ef_inverter,ef_battery);
        otherwise
            eBateria=descarga(SOC,ind_auto,e_gen,e_demand,ef_inverter);
            eBateria(eBateria>SOCMin)=SOCMin; %Porque la descarga es negativa
    end
 end