function [areaL,inverter,panel,turbina,battery,lco,clima,potencia_requerida,diesel]=cargaExcel()
fullYear=true; 
fullC=["25","26"];
if fullYear
    cantidad_datos=table2array(readtable("Entradas.xlsx","sheet","RecursoRenovable","range","A:A")); 
    auto=length(cantidad_datos(~isnan(cantidad_datos)));
    fullC=[string(auto+1),string(auto+2)];
end
ro=@(temperatura,h)(354.049./temperatura.*exp(-0.034.*h./temperatura));
rangos=["B2:B4";"B7:B9";"B12:B22";"B25:B32";"I3:I6";"B35:B35";"M3:M6";"M8:M10";"M12:M14";"M16:M21"]; 
nRangos=length(rangos); 
exportar=struct("a",{});
%%
for i=1:nRangos
    opts=spreadsheetImportOptions("NumVariables",1);
    opts.Sheet="input";
    opts.VariableTypes="double";
    opts.DataRange=rangos(i);
    reclamar=readtable("Entradas.xlsx",opts);
    exportar(i).a=reclamar; 
    reclamar_matrix=table2array(reclamar);
    switch i
        case 1
            areaL=reclamar_matrix(1);
            inverter.eficiencia=reclamar_matrix(2);
            clima.altura=reclamar_matrix(3);
        case 2
            panel.potencia=reclamar_matrix(1);
            panel.area=reclamar_matrix(3);
            opts.Sheet="Eficiencia del modulo FV";
            opts.DataRange=strcat("A2:A",fullC(1)); 
            panel.eficienciaOriginal=reclamar_matrix(2);
            panel.eficiencia=table2array(readtable("Entradas.xlsx",opts));
        case 3
            turbina.potencia=reclamar_matrix(1);
            turbina.eficiencia=reclamar_matrix(2);
            turbina.velocidadNominal=reclamar_matrix(3);
            turbina.velocidadArranque=reclamar_matrix(4);
            turbina.velocidadMaxima=reclamar_matrix(5);
            turbina.areaBarrido=reclamar_matrix(6);
            clima.densidadAire=reclamar_matrix(7);
            clima.densidad=reclamar_matrix(7);
            turbina.alturaReferencia=reclamar_matrix(8);
            turbina.alturaUsada=reclamar_matrix(9);
            turbina.alpha=reclamar_matrix(10);
            turbina.areaOcupada=reclamar_matrix(11);
        case 4
            battery.eficiencia=reclamar_matrix(1);
            battery.autoDescarga=reclamar_matrix(2);
            battery.SOCMax=reclamar_matrix(3);
            battery.SOCMin=reclamar_matrix(4);
            battery.amperioHora=reclamar_matrix(5);
            battery.voltaje=reclamar_matrix(6);
            battery.profDescarga=reclamar_matrix(7);
            battery.diasAuto=reclamar_matrix(8);
        case 5
            lco.sol=reclamar_matrix(1);
            lco.viento=reclamar_matrix(2);
            lco.bat=reclamar_matrix(3);
            lco.diesel=reclamar_matrix(4);
        case 6
            diesel.potencia=reclamar_matrix(1);
        case 7
            %Para costo nivelado de la energía
            panel.costo=reclamar_matrix(1);
            panel.coym=reclamar_matrix(2);
            inverter.costo=reclamar_matrix(3);
            panel.vidaUtil=reclamar_matrix(4);
            
        case 8
            %Para eolico
            turbina.costo=reclamar_matrix(1);
            turbina.coym=reclamar_matrix(2);
            turbina.vidaUtil=reclamar_matrix(3);
            
        case 9
            %Para Baterias
            battery.costo=reclamar_matrix(1);
            battery.coym=reclamar_matrix(2);
            battery.vidaUtil=reclamar_matrix(3);
        case 10
            %Para generador diesel
            diesel.costo=reclamar_matrix(1);
            diesel.coym=reclamar_matrix(2);
            diesel.vidaUtil=reclamar_matrix(3);
            diesel.hr=reclamar_matrix(4);
            diesel.precioCombustible=reclamar_matrix(5);
            diesel.vidaUtil=reclamar_matrix(5);
            diesel.valorSalvamento=reclamar_matrix(6);
    end
    
end

opts=spreadsheetImportOptions("NumVariables",3);
opts.Sheet="RecursoRenovable";
opts.DataRange=strcat("A3:C",fullC(2)); 
opts.VariableTypes=["double","double","double"];
ambiente=readtable("Entradas.xlsx",opts);
potencia_requerida=ambiente.Var1;
clima.irradiancia=ambiente.Var2;
clima.velViento=ambiente.Var3;
opts=spreadsheetImportOptions("NumVariables",1);
opts.Sheet="Eficiencia del modulo FV";
opts.DataRange= strcat("C2:C",fullC(1));
opts.VariableTypes="double";
temperatura=table2array(readtable("Entradas.xlsx",opts))+273.15;
clima.densidadAire=ro(temperatura,clima.altura);

end