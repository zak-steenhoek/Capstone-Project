% MAC Function
function [X,Y,AC,S,Macloc,quarterSweep,halfSweep] = MacCode(Root_Chord, Tip_Chord, Half_Span, Sweep_Angle)
    yloc_sweep = Half_Span*tan(Sweep_Angle);

    y1Root = Tip_Chord;
    y2Root = - Root_Chord - Tip_Chord;

    y1Tip = -yloc_sweep - Root_Chord - Tip_Chord;
    y2Tip = -yloc_sweep + Root_Chord;

    x = [0 Half_Span];
    y = [y1Root-y2Root y1Tip-y2Tip];
    x_val = x(1) - y(1) * (x(2) - x(1)) / (y(2) - y(1));
    y_val = y1Root + (y1Tip - y1Root) * (x_val - 0) / (Half_Span - 0);
    
    X = [0 Half_Span Half_Span 0 0];
    Y = [0 -yloc_sweep -yloc_sweep-Tip_Chord -Root_Chord 0];
    
    leadingEdgeMAC = X(1) + x_val * (Y(2) - Y(1)) / (X(2) - X(1));
    trailingEdgeMAC = X(4) + x_val * (Y(3) - Y(4)) / (X(3) - X(4)) - Root_Chord;
    Macloc = [leadingEdgeMAC trailingEdgeMAC];
    
    AC = [x_val leadingEdgeMAC - (leadingEdgeMAC-trailingEdgeMAC)*0.25];
    S = trapz([0 Half_Span Half_Span 0], [0 -yloc_sweep -yloc_sweep-Tip_Chord -Root_Chord]);


    quarterSweepPos = [-Root_Chord*0.25 -yloc_sweep-(Tip_Chord*0.25)];
    quarterSweep = atan(( quarterSweepPos(1) - quarterSweepPos(2) ) / Half_Span) * 180/pi;
    halfsweepPos = [-Root_Chord*0.5 -yloc_sweep-(Tip_Chord*0.5)];
    halfSweep = atan(( halfsweepPos(1) - halfsweepPos(2) ) / Half_Span) * 180/pi;

    figure()
    plot(X,Y); hold on;
    grid on; axis equal;
    plot([0 Half_Span],[y1Root y1Tip],'k--');
    plot([0 Half_Span],[y2Root y2Tip],'k--');
    plot([0 Half_Span],[quarterSweepPos(1) quarterSweepPos(2)],'r--');
    plot([0 Half_Span],[halfsweepPos(1) halfsweepPos(2)],'b--');
    scatter(x_val,y_val,'ko')
end