      FUNCTION fx(kd,m,vx,vy)
      real         fx
      real         kd,m,vx,vy
      fx = -kd*sqrt(vx*vx+vy*vy)*(vx)/m
      return
      end

