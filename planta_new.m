function [panel,turbina,battery,diesel_gen,lco,potencia]=planta_new(clima,panel,turbina,inverter,battery,lco,potenciaRequeria,hora)
%%
%Ecuaciones usadas
%@(d_aire,area_barrido,eficiencia_turbina,vel_viento,numero_turbina)(1/2.*d_aire.*area_barrido.*eficiencia_turbina.*(vel_viento.^3).*numero_turbina.*10^-3 );
%Potencia producida por la turbina 
%     (((10^-3).*eficiencia_turbina.*d_aire.*area_barrido.*(vel_viento.^3).*numero_turbina)./(2));

p_Turbina=@(eficiencia_turbina,d_aire,area_barrido,vel_viento,numero_turbina)...
    (1/2.*d_aire.*area_barrido.*eficiencia_turbina.*(vel_viento.^3).*numero_turbina.*10^-3 );

%Potencia producida por los paneles

%@(eficiencia_panel,area,irradiancia,numero_panel)(irradiancia.*eficiencia_panel.*area.*numero_panel.*10^-3);
 p_Panel=@(eficiencia_panel,area,irradiancia,numero_panel)...
     (irradiancia.*eficiencia_panel.*area.*numero_panel.*10^-3);

%Potencia de salida Diesel
p_Diesel=@(potencia_requerida,pot_renovable)(potencia_requerida-pot_renovable);
%Costo nivelado de la energia
lcoe=@(lcoSol,lcoViento,lcoBat,lcoDiesel,energy_sol,energy_wind,energy_bat,energy_diesel)((lcoSol+lcoViento+lcoBat+lcoDiesel)./eGen);
lcoe_formula=@(lcoe_used,energy_used)(lcoe_used.*energy_used);
%%
%panel_run=panel.cantidad(panel.cantidad>0);
panel_run=panel.cantidad;
%turbina_run=turbina.cantidad(turbina.cantidad>0);
turbina_run=turbina.cantidad;
%Generacion de potencia
pv_gen=p_Panel(panel.eficiencia(hora),panel.area,clima.irradiancia(hora),panel_run);
tur_gen=p_Turbina(turbina.eficiencia,clima.densidadAire(hora),clima.velViento(hora),turbina_run,turbina.areaBarrido);
%Energia renovable total
renovable_gen=sum([pv_gen,tur_gen],2);
%Sistemas que cargan baterias
posGG=logical(renovable_gen>=potenciaRequeria(hora));
posUG=~posGG;
%%
%Carga
if isempty(battery.SOCi(posGG))==0
    posMM_1=logical(battery.SOCi<battery.SOCMax);
    posMM=logical(posMM_1.*posGG);
    if isempty(battery.SOCi(posMM))==0
        %Carga la bateria
        battery.SOCi(posMM)=bateria(inverter.eficiencia,battery.eficiencia,...
            battery.SOCi(posMM),battery.SOCMax,battery.SOCMin,battery.autoDescarga,renovable_gen(posMM),...
            potenciaRequeria(hora),"carga");
    end
end
%%
%Descarga
if isempty(battery.SOCi(posUG))==0
    posMm_1=logical(battery.SOCi>battery.SOCMin);
    posMm=logical(posMm_1.*posUG);
    if isempty(battery.SOCi(posMm))==0
        antes=battery.SOCi(posMm);
        battery.SOCi(posMm)=bateria(inverter.eficiencia,battery.eficiencia,battery.SOCi(posMm),battery.SOCMax,...
            battery.SOCMin,battery.autoDescarga,renovable_gen(posMm),potenciaRequeria(hora),"descarga");
        %La perdida de energia
         battery.SOCL(posMm,hora)=antes-battery.SOCi(posMm);
         %La energia que aporta a la energia renovable
        renovable_gen(posMm,:)=renovable_gen(posMm,:)+battery.SOCi(posMm);
    end
end
%%
diesel_gen=p_Diesel(potenciaRequeria(hora),renovable_gen);
diesel_gen(diesel_gen<0)=diesel_gen(diesel_gen<0).*0;
energia_generada=renovable_gen+diesel_gen;
%lco.total=lcoe(lco.sol,lco.viento,lco.bat,lco.diesel,energia_generada);
potencia.panel=pv_gen;
potencia.turbina=tur_gen;
potencia.diesel=diesel_gen;
potencias=[pv_gen,tur_gen,battery.SOCi,diesel_gen];
lcoe_using=[lco.sol,lco.viento,lco.bat,lco.diesel];
[col_p,~]=size(pv_gen);
lco_temporal=zeros(col_p,4);
for i=1:4
    IM=lcoe_formula(lcoe_using(i),potencias(:,i));
    IM_replace01=isnan(IM);
    IM_replace02=isinf(IM);
    IM(logical(IM_replace01))=0;
    IM(logical(IM_replace02))=0;
    lco_temporal(:,i)=IM;
end
lco.total=sum(lco_temporal,2)./energia_generada;
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
    end
    eBateria(eBateria>SOCMax)=SOCMax;
    eBateria(eBateria<SOCMin)=SOCMin;
 end