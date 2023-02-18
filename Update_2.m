function Phi = Update_2(phi,v,dt,dx,dy)
    opt = (v >= 0);
    Phi_1 = phi;
    Phi_2 = phi;
    Phi_1(3:end-1,1:end) = phi(3:end-1,1:end) - dt*v(3:end-1,1:end).*( (phi(4:end,1:end) - phi(2:end-2,1:end))/(2*dy) ...
        - (phi(4:end,1:end) - 3*phi(3:end-1,1:end) + 3*phi(2:end-2,1:end) - phi(1:end-3,1:end))/(6*dy) );
    Phi_2(2:end-2,1:end) = phi(2:end-2,1:end) - dt*v(2:end-2,1:end).*( (phi(3:end-1,1:end) - phi(1:end-3,1:end))/(2*dy) ...
        - (phi(4:end,1:end) - 3*phi(3:end-1,1:end) + 3*phi(2:end-2,1:end) - phi(1:end-3,1:end))/(6*dy) );
    
    Phi = Phi_1.*opt + Phi_2.*(1 - opt);
    Phi(1,1:end) = Phi(3,1:end);
    Phi(end,1:end) = Phi(end-2,1:end);
    %grids adjacent to the lateral boundaries
    Phi_1(2,1:end) = phi(2,1:end) - dt*v(2,1:end).*(phi(2,1:end) - phi(1,1:end))/dy;
    Phi_2(2,1:end) = phi(2,1:end) - dt*v(2,1:end).*(phi(3,1:end) - phi(2,1:end))/dy;
    Phi_1(end-1,1:end) = phi(end-1,1:end) - dt*v(end-1,1:end).*(phi(end-1,1:end) - phi(end-2,1:end))/dy;
    Phi_2(end-1,1:end) = phi(end-1,1:end) - dt*v(end-1,1:end).*(phi(end,1:end) - phi(end-1,1:end))/dy;
    
    Phi = Phi_1.*opt + Phi_2.*(1 - opt);
end

