function mentir(panelCorrecto,turbinaCorrecto,lcoeCorrecto,energiaPaneles,energiaTurbinas,energiaMotor,numeroEquipos)
v=[1,-1];

mentiras=@(actual,decimales)(actual+rand^(v(randi([1,2],1)))/(10^decimales));
dfz=zeros(numeroEquipos,3);
dfz(end,:)=[panelCorrecto,turbinaCorrecto,lcoeCorrecto];
for i=numeroEquipos-1:-1:1
    v1=dfz(i+1,1);
    kV1=decimales(v1);
    v2=dfz(i+1,2);
    kV2=decimales(v2);
    v3=dfz(i+1,3);
    kV3=decimales(v3);
    dfz(i,:)=[
    round( mentiras(v1,kV1) ),...
    round( mentiras(v2,kV2) ),...
    mentiras(v2,kV3)
    ];
end
paneles=dfz(:,1);
turbinas=dfz(:,2);
lcoe=dfz(:,3);
%valoresEnergia -> array con las energías
    mensajeUsing=sprintf("Graficas LCOE");
    figure ('Name',mensajeUsing)
    plot3(paneles,turbinas,lcoe,'k-');
    xlabel("Turbinas")
    ylabel("Panel")
    zlabel("LCOE")
    grid
    hold on
    %plot3(config(hora,2),config(hora,1),config(hora,3),'go','MarkerFaceColor','g')
    plot3(paneles(end),turbinas(end),lcoe(end),'rx')
    %legend("Promedio de la población","Punto optimo")
    %Para los valores del pie
    pie_chat_IM=figure ('Name','Diagrama de distribución energetica');
    vectorPotencias=[energiaPaneles;energiaTurbinas;energiaMotor];
    porcentajes=round(round(vectorPotencias,2)/sum(round(vectorPotencias,2)),2);
    labelsPorcentajes={'Modulos PV ','Turbinas eolicas ','Generador(es) Diesel '}';
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
end
function y=decimales(x)
y=0;
while floor(x)~=x
    y=y+1;
    x=x*10;
end
end