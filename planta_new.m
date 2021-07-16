function [panel,turbina,battery,diesel_gen,lco,potencia]=planta_new(clima,panel,turbina,inverter,battery,lco,potenciaRequeria,hora)
%%
p_Turbina=@(eficiencia_turbina,d_aire,area_barrido,vel_viento,numero_turbina)...
    (1/2.*d_aire.*area_barrido.*eficiencia_turbina.*(vel_viento.^3).*numero_turbina.*10^-3 );

p_Panel=@(eficiencia_panel,area,irradiancia,numero_panel)...
     (irradiancia.*eficiencia_panel.*area.*numero_panel.*10^-3);

p_Diesel=@(potencia_requerida,pot_renovable)(potencia_requerida-pot_renovable);

inversion=@(numeroEquipos,costoEquipo1,costoEquipo2)(numeroEquipos.*(costoEquipo1+costoEquipo2));
aci=@(inversion,tasaInteres,vidaUtil)( (inversion.*(1+tasaInteres).^vidaUtil)./( (1+tasaInteres).^vidaUtil -1 ) );
acCom=@(energyDiesel,efGenerator,hr,priceDiesel)(energyDiesel./efGenerator.*hr.*priceDiesel);
lcoe=@(aciV,acom,energiaGenerada,acOther)((aciV+acom+acOther)./energiaGenerada);

%lcoe=@(lcoSol,lcoViento,lcoBat,lcoDiesel,energy_sol,energy_wind,energy_bat,energy_diesel)((lcoSol+lcoViento+lcoBat+lcoDiesel)./eGen);
lcoe_formula=@(lcoe_used,energy_used)(lcoe_used.*energy_used);
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
    end
end
%%
diesel_gen=p_Diesel(potenciaRequeria(hora),renovable_gen);
diesel_gen(diesel_gen<0)=diesel_gen(diesel_gen<0).*0;
energia_generada=renovable_gen+diesel_gen;

%Para energía total igual a cero
energy=( pv_gen+tur_gen+diesel_gen+battery.SOCi-battery.SOCL(:,hora)~=potenciaRequeria(hora) );

panel.cantidad(energy)=NaN;
panel.cantidad=fillmissing(panel.cantidad,'previous');
pv_gen(energy)=NaN;
pv_gen=fillmissing(pv_gen,'previous');
panel.inversion=inversion(panel.cantidad,panel.costo,inverter.costo);
panel.aci=aci(panel.inversion,panel.tasaInteres,panel.vidaUtil);
panel.lcoe=lcoe(panel.aci,0.015,pv_gen,0);

turbina.cantidad(energy)=NaN;
turbina.cantidad=fillmissing(turbina.cantidad,'previous');
tur_gen(energy)=NaN;
tur_gen=fillmissing(tur_gen,'previous');
turbina.inversion=inversion(turbina.cantidad,turbina.costo,0);
turbina.aci=aci(turbina.inversion,turbina.tasaInteres,turbina.vidaUtil);
turbina.lcoe=lcoe(turbina.aci,0.015,tur_gen,0);

diesel_gen(energy)=NaN;
diesel_gen=fillmissing(diesel_gen,'previous');
diesel.inversion=inversion(diesel_gen./diesel.potencia,diesel.costo,0);
diesel.acCom=(diesel_gen,diesel.eficiencia,diesel.hr,diesel.precioCombustible);
diesel.aci=aci(diesel.inversion,diesel.tasaInteres,diesel.vidaUtil);
diesel.lcoe=lcoe(diesel.aci,0.02,diesel_gen,diesel.acCom);

battery.SOCi(energy)=NaN;
battery.SOCi=fillmissing(battery.SOCi,'previous');
battery.inversion=inversion(battery.SOCi./battery.SOCMax,battery.costo,0);
battery.aci=aci(battery.inversion,battery.tasaInteres,battery.vidaUtil);
battery.lcoe=lcoe(battery.aci,0.02,battery.SOCi,0);
battery.SOCL(energy,hora)=NaN;
battery.SOCL=fillmissing(battery.SOCL(:,hora),'previous');

%Fin para energía total igual a cero

potencia.panel=pv_gen;
potencia.turbina=tur_gen;
potencia.diesel=diesel_gen;
potencias=[pv_gen,tur_gen,battery.SOCi,diesel_gen];

%lcoe_using=[lco.sol,lco.viento,lco.bat,lco.diesel];

%Para calcular de nuevo los LCOE
lco_temporal=[panel.lcoe,turbina.lcoe,diesel.lcoe,battery.lcoe];
lco.total=sum(lco_temporal,2)./energia_generada;
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