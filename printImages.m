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
    figure ("Name","Grafico Criterio de Parada") %Por ahora dejemos esto asÃ­
    plot(best_lcoe,'b-o')
    hold on
    plot(memoria_lcoe,'ro','MarkerFaceColor','r')
    grid
    xlabel("Iteracion")
    ylabel("LCOE")
    legend('Mejor LCOE','LCOE promedio')
end