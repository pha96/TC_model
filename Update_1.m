function Phi = Update_1(phi,u,dt,dx,dy)
    opt = (u >= 0);
    Phi_1 = phi;
    Phi_2 = phi;
    Phi_1(1:end,3:end-1) = phi(1:end,3:end-1) - dt*u(1:end,3:end-1).*( (phi(1:end,4:end) - phi(1:end,2:end-2))/(2*dx) ...
        - (phi(1:end,4:end) - 3*phi(1:end,3:end-1) + 3*phi(1:end,2:end-2) - phi(1:end,1:end-3))/(6*dx) );
    Phi_2(1:end,2:end-2) = phi(1:end,2:end-2) - dt*u(1:end,2:end-2).*( (phi(1:end,3:end-1) - phi(1:end,1:end-3))/(2*dx) ...
        - (phi(1:end,4:end) - 3*phi(1:end,3:end-1) + 3*phi(1:end,2:end-2) - phi(1:end,1:end-3))/(6*dx) );
    
    Phi = Phi_1.*opt + Phi_2.*(1 - opt);
    Phi(1:end,1) = Phi(1:end,3);
    Phi(1:end,end) = Phi(1:end,end-2);
    %grids adjacent to the lateral boundaries
    Phi_1(1:end,2) = phi(1:end,2) - dt*u(1:end,2).*(phi(1:end,2) - phi(1:end,1))/dx;
    Phi_2(1:end,2) = phi(1:end,2) - dt*u(1:end,2).*(phi(1:end,3) - phi(1:end,2))/dx;
    Phi_1(1:end,end-1) = phi(1:end,end-1) - dt*u(1:end,end-1).*(phi(1:end,end-1) - phi(1:end,end-2))/dx;
    Phi_2(1:end,end-1) = phi(1:end,end-1) - dt*u(1:end,end-1).*(phi(1:end,end) - phi(1:end,end-1))/dx;
    
    Phi = Phi_1.*opt + Phi_2.*(1 - opt);        
end

