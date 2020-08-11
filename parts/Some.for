      Program someparts
      real         g,dt,t0,tn,t[allocatable](:),ro
      real         part0[allocatable](:,:),vwx,vwy
      real         vx[allocatable](:,:),vy[allocatable](:,:)
      real         x[allocatable](:,:),y[allocatable](:,:)
      real         dist[allocatable](:,:,:)
      integer      nt,i,j,np
      
      g = 9.81
      np = 2
      ro = 1.17
      t0 = 0

      tn = 1.0
      nt = 10000
      dt = (tn-t0)/nt
      vwx = 0.0
      vwy = -2.9
      allocate(part0(np,9),vx(np,nt),vy(np,nt),x(np,nt),y(np,nt),t(nt))
      allocate(dist(nt,np,np))
      
      do j = 1,nt
          do i = 1,np
              do k = 1,np
                  dist(j,i,k) = 0.0
              end do
          end do
      end do
      
      !initial data
      part0(1,1) = 4.0/3.0*3.1415926*(75.0/1000000)**3.0*2550.0
      part0(1,2) = 0.47
      part0(1,3) = 150.0/1000000.0
      part0(1,4) = 0.05
      part0(1,5) = 10.0
      part0(1,6) = 1.0
      part0(1,7) = 2.0*3.1415926/3.0
     
      part0(2,1) = 4.0/3.0*3.1415926*(75.0/1000000)**3.0*2550
      part0(2,2) = 0.47
      part0(2,3) = 150.0/1000000.0
      part0(2,4) = 0.05
      part0(2,5) = 10.0
      part0(2,6) = 1.0
      part0(2,7) = 3.1415926/6.0
      
      t(1) = t0
      do i = 1,np
          part0(i,9) = part0(i,2)*ro*3.1415926*part0(i,3)*part0(i,3)/8.0
      end do
      
      do i=1,np
          vx(i,1) = part0(i,6)*cos(part0(i,7))
          vy(i,1) = part0(i,6)*sin(part0(i,7))
          x(i,1) = part0(i,4)
          y(i,1) = part0(i,5)
      end do
      
      do j = 2,nt
          t(j) = t(j-1)+dt
          do i = 1,np   
              call RungK(g,vwx,vwy,part0(i,1:9),vx(i,j-1),vy(i,j-1),
     *                  x(i,j-1),y(i,j-1),
     *                  vx(i,j),vy(i,j),x(i,j),y(i,j),dt)
              if (y(i,j).LE.(part0(i,3)/2.0)) then
                  vy(i,j) = -vy(i,j)
              end if
              if (x(i,j).LE.(part0(i,3)/2.0).or.(x(i,j).GE.(-part0(i,3)
     *            /2.0 + 0.1))) then
                  vx(i,j) = -vx(i,j)
              end if
          end do
          
          do i = 1,np
          do k = i+1,np
                  dist(j,i,k) = sqrt(((y(1,j)-y(2,j))**2)
     *                     +(x(1,j)-x(2,j))**2)
                  !write(*,*)dist(j,i,k)
                  if (dist(j,i,k) .le. part0(i,3)/2.0 + part0(k,3)/2.0) 
     *            then
                      !write(*,*)'collision'
                      !write(*,*)t(j)
!                      vx(i,j)=vx(i,j)-part0(k,1)/(part0(k,1)
!     *                        +part0(i,1))*(vx(i,j)-vx(k,j))*1.0
!                      vy(i,j)=vy(i,j)-part0(k,1)/(part0(k,1)
!     *                        +part0(i,1))*(vy(i,j)-vy(k,j))*1.0
!                      vx(k,j)=vx(k,j)-part0(i,1)/(part0(i,1)
!     *                        +part0(k,1))*(vx(k,j)-vx(i,j))*1.0
!                      vy(k,j)=vy(k,j)-part0(i,1)/(part0(i,1)
!     *                        +part0(k,1))*(vy(k,j)-vy(i,j))*1.0
                 end if
          end do
          end do
      end do
      
      open(5,file='xy.dat')
      open(10,file='vxvy.dat')
      do j = 1,nt
      write(5,*)t(j),x(1,j),y(1,j),x(2,j),y(2,j)
      write(10,*)t(j),vx(1,j),vy(1,j),vx(2,j),vy(2,j)
      end do
      
      pause
      end program