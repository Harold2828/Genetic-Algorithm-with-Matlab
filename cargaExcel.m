function [areaL,inverter,panel,turbina,battery,lco,clima,potencia_requerida]=cargaExcel()
%LMPC
fullYear=true; %<---
fullC=["25","26"];
if fullYear
    auto=length(table2array(readtable("Entradas.xlsx","sheet","Eficiencia del modulo FV","range","A:A")));
    fullC=[string(auto+1),string(auto+2)];
end
ro=@(temperatura,h)(354.049./temperatura.*exp(-0.034.*h./temperatura));
rangos=["B2:B4";"B7:B9";"B12:B19";"B22:B27";"I3:I6"]; 
nRangos=length(rangos); 
exportar=struct("a",{});
%%
%sheet 1 -> Input
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
            %Datos generales
            areaL=reclamar_matrix(1);
            inverter.eficiencia=reclamar_matrix(2)./100;
            clima.altura=reclamar_matrix(3);
        case 2
            %Modelo fotovoltaico
            panel.potencia=reclamar_matrix(1);
            panel.area=reclamar_matrix(3);
            opts.Sheet="Eficiencia del modulo FV";
            opts.DataRange=strcat("A2:A",fullC(1)); 
            panel.eficiencia=table2array(readtable("Entradas.xlsx",opts));
        case 3
            %Modelo eolico
            turbina.potencia=reclamar_matrix(1);
            turbina.eficiencia=reclamar_matrix(2)./100;
            turbina.areaBarrido=reclamar_matrix(3);
            clima.densidadAire=reclamar_matrix(4);
            turbina.alturaReferencia=reclamar_matrix(5);
            turbina.alturaUsada=reclamar_matrix(6);
            turbina.alpha=reclamar_matrix(7);
            turbina.areaOcupada=reclamar_matrix(8);
        case 4
            %Modelo baterias
            battery.eficiencia=reclamar_matrix(1)./100;
            battery.autoDescarga=reclamar_matrix(2);
            battery.SOCMax=reclamar_matrix(3);
            battery.SOCMin=reclamar_matrix(4);
            battery.profDescarga=reclamar_matrix(5);
            battery.diasAuto=reclamar_matrix(6);
        case 5
            %Costo nivelado de la energía
            lco.sol=reclamar_matrix(1);%LMPC
            lco.viento=reclamar_matrix(2);
            lco.bat=reclamar_matrix(3);
            lco.diesel=reclamar_matrix(4);
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
%Para lo de la densidad
opts=spreadsheetImportOptions("NumVariables",1);
opts.Sheet="Eficiencia del modulo FV";
opts.DataRange= strcat("C2:C",fullC(1));
opts.VariableTypes="double";
temperatura=table2array(readtable("Entradas.xlsx",opts))+273.15;
clima.densidadAire=ro(temperatura,clima.altura);

end