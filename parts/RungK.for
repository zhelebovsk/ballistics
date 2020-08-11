      Subroutine RungK(g,vwx,vwy,part0,vx0,vy0,x0,y0,vx,vy,x,y,dt)
      real    part0(9),vx0,vy0,vx,vy,x0,y0,x,y
      real    g,vwx,vwy,m,kd,dt
      real    kx(4),ky(4)
      
      kd = part0(9)
      m = part0(1)
      
      kx(1) = fx(kd,m,vx0+vwx,vy0+vwy)
      ky(1) = fy(kd,m,vx0+vwx,vy0+vwy,g)
      kx(2) = fx(kd,m,(vx0+vwx)+kx(1)*dt/2.0,(vy0+vwy)+ky(1)*dt/2.0)
      ky(2) = fy(kd,m,(vx0+vwx)+kx(1)*dt/2.0,(vy0+vwy)+ky(1)*dt/2.0,g)
      kx(3) = fx(kd,m,(vx0+vwx)+kx(2)*dt/2.0,(vy0+vwy)+ky(2)*dt/2.0)
      ky(3) = fy(kd,m,(vx0+vwx)+kx(2)*dt/2.0,(vy0+vwy)+ky(2)*dt/2.0,g)
      kx(4) = fx(kd,m,(vx0+vwx)+kx(3)*dt,(vy0+vwy)+ky(3)*dt)
      ky(4) = fy(kd,m,(vx0+vwx)+kx(3)*dt,(vy0+vwy)+ky(3)*dt,g)
      if (x0 > 0.0) then 
          vx = vx0 + dt*(kx(1) + 2.0*kx(2) + 2.0*kx(3) + kx(4))/6.0
      else
          vx = vx0 - dt*(kx(1) + 2.0*kx(2) + 2.0*kx(3) + kx(4))/6.0
      end if
      if (y0 > 0.0) then 
          vy = vy0 + dt*(ky(1) + 2.0*ky(2) + 2.0*ky(3) + ky(4))/6.0
      else
          vy = vy0 - dt*(ky(1) + 2.0*ky(2) + 2.0*ky(3) + ky(4))/6.0
      end if    
      x = x0 + vx0*dt
      y = y0 + vy0*dt
      
      return
      end