function MI_AL_plot(i,ObjVal,ConViol,maxiter)
    figure(1);
    if i == 1
        set(gcf,'position',[10 30 400 300])
        hold on;
        grid on;
        title('Objective Function Value')
        xlabel('Iterations')
    end
    plot(i,ObjVal,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);
    
    figure(2);
    if i == 1    
        set(gcf,'position',[420 30 400 300])
        hold on;
        grid on;
        title('Equality Constraint Violation')
        xlabel('Iterations')
        axis([1 maxiter 0 1])  
    end
    plot(i,ConViol,'-mo','MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63]);

end