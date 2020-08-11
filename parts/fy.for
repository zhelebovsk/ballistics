      Function fy(kd,m,vx,vy,g)
      real         fy
      real         kd,m,vx,vy,g
      fy = -kd*sqrt(vx*vx+vy*vy)*(vy)/m - g
      return
      end

