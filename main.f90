      !=================================================================
      program main
         implicit none
         !gluons
         !1:x, 2:y, 3:z, 4:tau, 5:px, 6:py, 7:pz, 8:m, 9: not used
         real(8) :: g(9, 200000) ! gluons
         real(8) :: c(9, 200000) ! gluons
         real(8) :: cbar(9, 200000) ! gluons
         real(8) :: jpsi(9, 200000) ! gluons
         real(8) :: gid(200000) ! gluons

         real(8) :: mass_number
         real(8) :: radius, radius_N
         real(8) :: Temp, dz
         real(8) :: pi
         integer(8) :: n_gluon, n_pair, n_group, n_c
         integer(4) :: nseed
         integer(8) :: sig_id_gg, sig_id_cg ! 1:identical particle. 0: different particle
         integer(8) :: sig_update_gg, sig_update_cg
         integer(8) :: sig_cycl, sig_corr
         real(8) :: dt ! time step
         real(8) :: n_coll_gg, n_coll_cg, n_coll_cbarg
         real(8) :: n_coll_gg_th, n_coll_cg_th
         real(8) :: sigma_gg, sigma_cg, sigma_cc, sigma_jg
         real(8) :: t, t0, t_max
         real(8) :: lx, ly, lz
         real(8) :: x_up, y_up, z_up
         real(8) :: x_low, y_low, z_low
         integer(8) :: i_time
         real(8) :: m_c, m_psi ! mass of charm
         real(8) :: w_jpsi_pro, w_jpsi_sur
         integer(8) :: n_jpsi_pro, n_jpsi_sur
         real(8) :: volume_0 ! initial volume 
         real(8) :: fmGeVfactor

         fmGeVfactor = 0.197327

         mass_number = 208
         radius_N    = 1.24
         dz          = 0.5
         nseed       = 41  ! seed for random number
         sigma_gg    = 0.4 ! cross section of g+g->g+g in fm^2
         sigma_cg    = 0.4 ! cross section of c+g->c+g in fm^2
         sigma_cc    = 0.4 ! cross section of c+c->Jpsi+g in fm^2
         sigma_jg    = 0.4 ! cross section of Jpsi+g-> c+c in fm^2
         dt = 0.2
         t0 = 0.0
         t_max = 30.0
         sig_id_gg     = 1 ! identical for gg collision
         sig_update_gg = 3
         sig_id_cg     = 0 ! identical for gg collision
         sig_update_cg = 1
         sig_cycl = 1
         sig_corr = 0
         lx = 10.0
         ly = 10.0
         lz = 10.0
         x_low = -lx/2.0
         y_low = -ly/2.0
         z_low = -lz/2.0
         x_up  =  lx/2.0
         y_up  =  ly/2.0
         z_up  =  lz/2.0
         m_c   = 1.25
         m_psi = 3.1

         pi = 3.14159265359
         
         n_gluon = 1000
         Temp = 0.15
         n_gluon = int(40*lx*ly*lz*Temp**3/pi**2/fmGeVfactor**3)
!         Temp    = 0.3
         n_pair  = 1000
         n_group = 1

         !initialization
         call SRAND(nseed)
         radius = radius_N*mass_number**(1.0/3.0)
         print *, "radius=", radius, "fm"
         volume_0 = pi*radius*radius*dz
         !Temp = (n_gluon/volume_0*fmGeVfactor**3*pi**2/40.0)**(1.0/3.0)
         print *, "Temp=", Temp
         print *, "N_g=", N_gluon
         n_jpsi_sur=0
         n_jpsi_pro=0
         i_time = 0
         t      = t0
         n_c  = n_pair*n_group
         n_coll_gg_th=n_gluon*n_gluon*sigma_gg*dt/(2.0*lx*ly*lz)
         n_coll_cg_th=n_gluon*n_c*sigma_cg*dt/(lx*ly*lz)

         call sr_initialize_glon(g, n_gluon, Temp, dz, radius)
         !call sr_output_particle(g, n_gluon)

         call sr_initialize_c(c, cbar, gid, n_pair, n_group,  m_c,&
                              dz, radius, sig_corr)
         !Cascade

 1000    continue
!         call sr_jpsi_production(c, cbar, jpsi, gid, n_pair, n_group,&
!                                    n_jpsi_pro, n_jpsi_sur, m_psi, &
!                                    sigma_cc, dt, sig_corr, t)
!         call sr_jpsi_dissociation(jpsi, g, n_jpsi_pro, n_jpsi_sur,&
!                                      n_gluon, m_psi, m_c, sigma_jg, dt)
         call sr_collision_one_step_el(c, g, n_c, n_gluon,&
          sigma_cg, dt, sig_id_cg, n_coll_cg, sig_update_cg, t)
         call sr_collision_one_step_el(cbar, g, n_c, n_gluon,&
          sigma_cg, dt, sig_id_cg, n_coll_cbarg, sig_update_cg, t)
         call sr_collision_one_step_el(g, g, n_gluon, n_gluon,&
          sigma_gg, dt, sig_id_gg, n_coll_gg, sig_update_gg, t)

         call sr_jpsi_pro (jpsi, n_jpsi_pro, w_jpsi_pro)
         call sr_jpsi_sur (jpsi, n_jpsi_pro, w_jpsi_sur)

 901     format(8f15.8)
         write(*,901) t, n_coll_gg, n_coll_gg_th,&
                  n_coll_cg, n_coll_cbarg, n_coll_cg_th,&
                  w_jpsi_sur, w_jpsi_pro

         call sr_free_streaming(g, n_gluon, dt, x_up, x_low, &
            y_up, y_low, z_up, z_low, sig_cycl)
         call sr_free_streaming(c, n_c, dt, x_up, x_low, &
            y_up, y_low, z_up, z_low, sig_cycl)
         call sr_free_streaming(cbar, n_c, dt, x_up, x_low, &
            y_up, y_low, z_up, z_low, sig_cycl)
!         call sr_free_streaming(jpsi, n_jpsi_pro, dt, x_up, x_low, &
!            y_up, y_low, z_up, z_low, sig_cycl)

         i_time = i_time+1
         t = t0+i_time*dt         
         if(t .lt. t_max) goto 1000

         !call sr_output_particle(g, n_gluon)
         print *, "Done" 
      end program main
      !=================================================================


      !=================================================================
      ! initialize_glons
      !-----------------------------------
      subroutine sr_initialize_glon(g, n_gluon, Temp, dz, radius)
         implicit none
         real(8), intent(out) :: g(9, 200000)
         real(8), intent(in) :: Temp, dz, radius
         integer(8), intent(in) :: n_gluon

         integer(8) :: i
         real(8) :: x, y, z, px, py, pz, m
         real(8) :: rho, phi
         real(8) :: pi
         real(8) :: rand1, rand2, f, p
         real(8) :: cos_theta, sin_theta 
         pi = atan(1.0)*4.0
         m  = 0.0

         do i=1, n_gluon
 
            z = dz*(rand()-0.5)
            rho = sqrt(rand())*radius
            phi = rand()*2.0*pi
            x = rho*cos(phi)
            y = rho*sin(phi)
            
  200       rand1 = 15.0*rand()
            f = rand1**2*exp(-rand1)
            rand2 = 0.6*rand()

            if (rand2 .gt. f) goto 200

            p = Temp*rand1
            cos_theta = 2.0*rand()-1.0
            sin_theta = sqrt(1.0-cos_theta**2)
            phi       = 2.0*pi*rand() 

            px = p*sin_theta*cos(phi)
            py = p*sin_theta*sin(phi)
            pz = p*cos_theta

            g(1, i) = x
            g(2, i) = y
            g(3, i) = z
            g(4, i) = -0.1
            g(5, i) = px
            g(6, i) = py
            g(7, i) = pz 
            g(8, i) = m
            g(9, i) = 0.0
         enddo
      end subroutine


      !=================================================================
      ! inner product of 4-vectors A and B
      !-----------------------------------
      subroutine sr_dot ( A, B, res)
         implicit none
         real(8), intent(in)  :: A(4), B(4)
         real(8), intent(out) :: res
         res = A(1)*B(1)-A(2)*B(2)-A(3)*B(3)-A(4)*B(4)
      end subroutine sr_dot
      !=================================================================

      !=================================================================
      ! to check whether two particles collides
      ! time components x1(1) and x2(1) are assumed to be the same
      !-----------------------------------
      subroutine sr_coll_check (x1, x2, p1, p2, dt, sigma, sig_coll,&
         dt_coll)
         implicit none
         real(8), intent(in) :: x1(4), x2(4), p1(4), p2(4)
         real(8), intent(in) :: sigma, dt
         integer(8), intent(out) :: sig_coll
         real(8), intent(out) :: dt_coll

         real(8) :: dx(4)
         real(8) :: dxp1, dxp2, p12, dx_2
         real(8) :: dxp1_2, dxp2_2, p12_2, m1_2, m2_2, m12_2

         real(8) :: b2
         real(8) :: pi

         real(8) :: dis2_max

         pi  = 3.1415926535
         sig_coll = 0

         dx  = x2-x1
         dis2_max = sigma/pi + 4.0*dt*dt

         call sr_dot(dx, dx, dx_2)
         if(-dx_2 .gt. dis2_max) goto 100

         call sr_dot(dx, p1, dxp1)
         call sr_dot(dx, p2, dxp2)
         call sr_dot(p1, p2, p12)
         call sr_dot(p1, p1, m1_2)
         call sr_dot(p2, p2, m2_2)

         dxp1_2 = dxp1*dxp1
         dxp2_2 = dxp2*dxp2
         p12_2  = p12*p12
         m12_2  = m1_2*m2_2

         b2 = -dx_2-(dxp1_2*m2_2+dxp2_2*m1_2-2.0*dxp1*dxp2*p12)&
              /(p12_2-m12_2)
         dt_coll = (dxp2*(p12*p1(1)+m1_2*p2(1))-dxp1*(p12*p2(1)+&
                   m2_2*p1(1)))/(2.0*(p12_2-m12_2))

         if((b2*pi .lt. sigma) .and. (dt_coll .gt. 0.0) .and.&
         (dt_coll .lt. dt))then
            sig_coll = 1
         endif
  100    continue
      endsubroutine
            

      !=================================================================
      ! Test dot
      !-----------------------------------
      subroutine sr_test_dot ()
         real(8) :: x(4), y(4)
         real(8) :: res
         print *, "Hello, World!"
         x(1) = 1.0
         x(2) = 2.0
         x(3) = 3.0
         x(4) = 4.0
         y(1) = 5.0
         y(2) = 6.0
         y(3) = 7.0
         y(4) = 8.0
         call sr_dot(x,y,res)
         print *, "x.y=", res
      endsubroutine
      !=================================================================
     
      subroutine sr_output_particle( par, length)
         implicit none
         real(8)    :: par(9, 200000)
         integer(8) :: length

         integer(8) :: i
 2001    format(9f10.5)
         do i=1, length
            write(*, 2001) par(1,i), par(2,i), par(3,i), par(4,i),&
            par(5,i), par(6,i), par(7,i), par(8,i), par(9,i)
         enddo
      endsubroutine


      subroutine xLorentz(t1,x1,y1,z1, vx,vy,vz, t2,x2,y2,z2)
         implicit real*8 (a-h,o-z)
     
         beta2=vx*vx +vy*vy +vz*vz
         if (beta2.gt.1.d-9) then  !<----------gai
            gam=1./sqrt(1.-beta2)
        
            t2=gam*(t1 -vx*x1 -vy*y1 -vz*z1)
            x2=x1 -vx*gam*t1 +(gam-1.)*vx*vx/beta2*x1 &
               +(gam-1.)*vx*vy/beta2*y1 +(gam-1.)*vx*vz/beta2*z1
            y2=y1 -vy*gam*t1 +(gam-1.)*vy*vx/beta2*x1 &
               +(gam-1.)*vy*vy/beta2*y1 +(gam-1.)*vy*vz/beta2*z1
            z2=z1 -vz*gam*t1 +(gam-1.)*vz*vx/beta2*x1 &
               +(gam-1.)*vz*vy/beta2*y1 +(gam-1.)*vz*vz/beta2*z1
         else
            t2=t1
            x2=x1
            y2=y1
            z2=z1
         end if
     
         return
      end

      subroutine sr_collision_one_step_el(particle_1, particle_2, length_1, length_2,&
          sigma, dt, sig_id, n_coll, sig_update, time)

         implicit none

         !parameter variables=================================================
         ! array for the first species of particles
         ! The 8 components: 1-x, 2-y, 3-z, 4-t, 5-px, 6-py, 7-pz, 8-m
         real(8), intent(inout), dimension(9,200000) :: particle_1 
  
         ! array for the second specie of particles
         real(8), intent(inout), dimension(9,200000) :: particle_2
  
         real(8), intent(in) :: sigma ! total cross section
         real(8), intent(in) :: dt    ! time step
         integer(8), intent(in) :: length_1 ! number of particles in particle_1
         integer(8), intent(in) :: length_2 ! number of particles in particle_2
  
         ! signal for identical partical of specie 1 and 2.
         ! sig_id = 1 for identical particle, and 0 for different particle
         integer(8), intent(in) :: sig_id
  
         real(8), intent(out) :: n_coll ! collision rate
  
         ! for periodic boundary condition
         integer(8), intent(in) :: sig_update
         real(8), intent(in) :: time
         !====================================================================
  
         !temprary variables==================================================
         real(8) :: t  ! time 
  
         integer(8) :: i_min_1, i_min_2, i_max_1, i_max_2 !Loop boundary
         integer(8) :: i_1, i_2
         real(8) :: x_1, y_1, z_1, px_1, py_1, pz_1, m_1, E_1 ! particle_1
         real(8) :: x_2, y_2, z_2, px_2, py_2, pz_2, m_2, E_2 ! particle_2
  
         !_m for meeting (minimum distance in CMF)
  
         real(8) :: dis2, dis2_max, dis_min
         real(8) :: pi
         real(8) :: epsilon_v2
         real(8) :: s  ! square of the total energy in the center of mass frame
         real(8) :: ds ! ds= s-(m_1+m_2)^2, ds+4*m_1*m_2=s-(m_1-m_2)^2
         real(8) :: pr2 ! momentum square of one particle in the CMF
         real(8) :: pr  ! pr=sqrt(pr2)
  
         real(8) :: x1(4), x2(4), x3(4), x4(4)
         real(8) :: p1(4), p2(4), p3(4), p4(4), p(4)
         real(8) :: v1(4), v2(4), v3(4), v4(4), v(4)
  
         real(8) :: E_c_3, px_c_3, py_c_3, pz_c_3
         real(8) :: E_3,   px_3,   py_3,   pz_3
  
         real(8) :: dt_coll
         integer(8) :: sig_coll
         real(8) :: cos_theta, sin_theta, phi ! random angle after collision 
  
         ! 0: update neither, 1: only update particle_1, 
         ! 2: only update particle_2, 3: update both
         !integer(8) :: sig_update
  
         !====================================================================
         !initialization======================================================
         pi = atan(1.0)*4.0
         dis2_max = sigma/pi + 4.0*dt*dt ! An upper boundary to collide 
         dis_min = 1.d-9 ! very small distance is not allowed
         !sig_update = sig_update_0
  
         i_min_1 = 1
         if(sig_id .eq. 1) then
            i_max_1=length_1-1
         else
            i_max_1=length_1
         endif
         epsilon_v2 = 1.d-14 ! to avoid v>=c
         n_coll = 0.0
         t = 0.0
         !====================================================================
  
  
         !core code===========================================================
         do i_1 = i_min_1, i_max_1
 
            if(time .lt. particle_1(4,i_1)) goto 3001 
            x_1  = particle_1(1,i_1)
            y_1  = particle_1(2,i_1)
            z_1  = particle_1(3,i_1)
  
            px_1 = particle_1(5,i_1)
            py_1 = particle_1(6,i_1)
            pz_1 = particle_1(7,i_1)
            m_1  = particle_1(8,i_1)
  
            E_1  = sqrt(px_1**2+py_1**2+pz_1**2+m_1**2)
  
            p1(1)= E_1
            p1(2)= px_1
            p1(3)= py_1
            p1(4)= pz_1
            x1(1)= 0.0
            x1(2)= x_1
            x1(3)= y_1
            x1(4)= z_1
  
            v1 = p1/E_1
  
            if(sig_id .eq. 1) then
               i_min_2 = i_1+1
            else
               i_min_2 = 1
            endif
            i_max_2 = length_2 
  
            do i_2 = i_min_2, i_max_2
               
               if(time .lt. particle_2(4,i_2)) goto 3000 
               x_2  = particle_2(1,i_2)
               y_2  = particle_2(2,i_2)
               z_2  = particle_2(3,i_2)
  
               dis2 = (x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2)&
                     +(z_1-z_2)*(z_1-z_2)
               if(dis2 .gt. dis2_max) goto 3000
     
               px_2 = particle_2(5,i_2)
               py_2 = particle_2(6,i_2)
               pz_2 = particle_2(7,i_2)
               m_2  = particle_2(8,i_2)
  
               E_2 = sqrt(px_2**2+py_2**2+pz_2**2+m_2**2)
  
               p2(1)= E_2
               p2(2)= px_2
               p2(3)= py_2
               p2(4)= pz_2
               x2(1)= 0.0
               x2(2)= x_2
               x2(3)= y_2
               x2(4)= z_2
  
               call sr_coll_check (x1, x2, p1, p2, dt, sigma,&
                                   sig_coll, dt_coll)
               if(sig_coll .eq. 0) goto 3000
               p = p1+p2
               call sr_dot(p,p,s)
               ds = s-(m_1+m_2)**2 
!             ds = s-(m_3+m_4)**2 
!             if( ds .le. 0.0) goto 3000

               ! Collision happens
               n_coll = n_coll +1  
  
               pr2 = ds*(ds+4.0*m_1*m_2)/(4.0*s)
               pr  = sqrt(pr2)
  
               v = p/p(1) !velocity (1, vx, vy, vz) of center of mass
  
               cos_theta = 2.0*rand()-1.0
               sin_theta = sqrt(1.0-cos_theta*cos_theta)
               phi       = 2.0*pi*rand()
  
               pz_c_3 = pr*cos_theta
               px_c_3 = pr*sin_theta*cos(phi)
               py_c_3 = pr*sin_theta*sin(phi)
               E_c_3  = sqrt(pr2+m_1**2)
  
               call xLorentz(E_c_3, px_c_3, py_c_3, pz_c_3, -v(2), -v(3),&
                    -v(4), E_3, px_3, py_3, pz_3)
               p3(1) = E_3
               p3(2) = px_3
               p3(3) = py_3
               p3(4) = pz_3
               v3 = p3/E_3
  
               x3 = x1+(v1-v3)*dt_coll
  
               p4 = p-p3
               v4 = p4/p4(1)
               x4 = x2+(v2-v4)*dt_coll
  
               if ((sig_update .eq. 1) .or. (sig_update .eq. 3)) then 
  
                  particle_1(1, i_1) = x3(2)
                  particle_1(2, i_1) = x3(3)
                  particle_1(3, i_1) = x3(4)
                  particle_1(5, i_1) = p3(2)
                  particle_1(6, i_1) = p3(3)
                  particle_1(7, i_1) = p3(4)
          
               endif
  
               if ((sig_update .eq. 2) .or. (sig_update .eq. 3)) then 
  
                  particle_2(1, i_2) = x4(2)
                  particle_2(2, i_2) = x4(3)
                  particle_2(3, i_2) = x4(4)
                  particle_2(5, i_2) = p4(2)
                  particle_2(6, i_2) = p4(3)
                  particle_2(7, i_2) = p4(4)
       
               endif
  
               goto 3001
           
 3000          continue             

            enddo

 3001       continue

         enddo

      endsubroutine sr_collision_one_step_el

      subroutine sr_jpsi_production(c, cbar, jpsi, gid, n_pair, n_group,&
                                    n_jpsi_pro, n_jpsi_sur, m_psi, &
                                    sigma, dt, sig_corr, time)

         implicit none

         !parameter variables=================================================
         ! array for the first species of particles
         ! The 8 components: 1-x, 2-y, 3-z, 4-t, 5-px, 6-py, 7-pz, 8-m
         real(8), intent(inout), dimension(9,200000) :: c
         real(8), intent(inout), dimension(9,200000) :: cbar
         real(8), intent(inout), dimension(9,200000) :: jpsi
         real(8), intent(inout), dimension(200000) :: gid
  
         real(8), intent(in) :: sigma ! total cross section
         real(8), intent(in) :: dt    ! time step
         real(8), intent(in) :: m_psi
         integer(8), intent(in) :: n_pair, n_group ! number of particles in particle_1
         integer(8), intent(inout) :: n_jpsi_pro, n_jpsi_sur
         integer(8), intent(in) :: sig_corr
         integer(8), intent(in) :: time
  
         !====================================================================
  
         !temprary variables==================================================
         real(8) :: t  ! time 
         integer(8) :: n_c
  
         integer(8) :: i_1, i_2
         real(8) :: x_1, y_1, z_1, px_1, py_1, pz_1, m_1, E_1 ! particle_1
         real(8) :: x_2, y_2, z_2, px_2, py_2, pz_2, m_2, E_2 ! particle_2
  
         !_m for meeting (minimum distance in CMF)
  
         real(8) :: dis2, dis2_max, dis_min
         real(8) :: pi
         real(8) :: epsilon_v2
         real(8) :: s  ! square of the total energy in the center of mass frame
         real(8) :: ds ! ds= s-(m_1+m_2)^2, ds+4*m_1*m_2=s-(m_1-m_2)^2
         real(8) :: pr2 ! momentum square of one particle in the CMF
         real(8) :: pr  ! pr=sqrt(pr2)
  
         real(8) :: x1(4), x2(4), x3(4)
         real(8) :: p1(4), p2(4), p3(4), p(4)
         real(8) :: v1(4), v3(4), v(4)
  
         real(8) :: E_c_3, px_c_3, py_c_3, pz_c_3
         real(8) :: E_3,   px_3,   py_3,   pz_3
  
         real(8) :: dt_coll
         integer(8) :: sig_coll
         real(8) :: cos_theta, sin_theta, phi ! random angle after collision 

         real(8) :: m_g, weight
  
         ! 0: update neither, 1: only update particle_1, 
         ! 2: only update particle_2, 3: update both
         !integer(8) :: sig_update
  
         !====================================================================
         !initialization======================================================
         pi = atan(1.0)*4.0
         dis2_max = sigma/pi + 4.0*dt*dt ! An upper boundary to collide 
         dis_min = 1.d-9 ! very small distance is not allowed
         m_g = 0.0
         n_c = n_pair*n_group
         !sig_update = sig_update_0
  
         epsilon_v2 = 1.d-14 ! to avoid v>=c
         t = 0.0
         !====================================================================
  
  
         !core code===========================================================
         do i_1 = 1, n_c
  
            if(time .lt. c(4,i_1)) goto 3001            
            x_1  = c(1,i_1)
            y_1  = c(2,i_1)
            z_1  = c(3,i_1)
  
            px_1 = c(5,i_1)
            py_1 = c(6,i_1)
            pz_1 = c(7,i_1)
            m_1  = c(8,i_1)
  
            E_1  = sqrt(px_1**2+py_1**2+pz_1**2+m_1**2)
  
            p1(1)= E_1
            p1(2)= px_1
            p1(3)= py_1
            p1(4)= pz_1
            x1(1)= 0.0
            x1(2)= x_1
            x1(3)= y_1
            x1(4)= z_1
  
            v1 = p1/E_1
  
            do i_2 = 1, n_c

               if((sig_corr.eq.1).and.(gid(i_1).eq.gid(i_2))&
                  .and.(i_1.ne.i_2)) goto 3000
               if(time .lt. cbar(4,i_2)) goto 3000            

               x_2  = cbar(1,i_2)
               y_2  = cbar(2,i_2)
               z_2  = cbar(3,i_2)
  
               dis2 = (x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2)&
                     +(z_1-z_2)*(z_1-z_2)
               if(dis2 .gt. dis2_max) goto 3000
     
               px_2 = cbar(5,i_2)
               py_2 = cbar(6,i_2)
               pz_2 = cbar(7,i_2)
               m_2  = cbar(8,i_2)
  
               E_2 = sqrt(px_2**2+py_2**2+pz_2**2+m_2**2)
  
               p2(1)= E_2
               p2(2)= px_2
               p2(3)= py_2
               p2(4)= pz_2
               x2(1)= 0.0
               x2(2)= x_2
               x2(3)= y_2
               x2(4)= z_2
  
               call sr_coll_check (x1, x2, p1, p2, dt, sigma,&
                                   sig_coll, dt_coll)
               if(sig_coll .eq. 0) goto 3000
               p = p1+p2
               call sr_dot(p,p,s)
               ds = s-(m_psi+m_g)**2 
               if( ds .le. 0.0) goto 3000

               ! Collision happens
  
               pr2 = ds*(ds+4.0*m_psi*m_g)/(4.0*s)
               pr  = sqrt(pr2)
  
               v = p/p(1) !velocity (1, vx, vy, vz) of center of mass
  
               cos_theta = 2.0*rand()-1.0
               sin_theta = sqrt(1.0-cos_theta*cos_theta)
               phi       = 2.0*pi*rand()
  
               pz_c_3 = pr*cos_theta
               px_c_3 = pr*sin_theta*cos(phi)
               py_c_3 = pr*sin_theta*sin(phi)
               E_c_3  = sqrt(pr2+m_psi**2)
  
               call xLorentz(E_c_3, px_c_3, py_c_3, pz_c_3, -v(2), -v(3),&
                    -v(4), E_3, px_3, py_3, pz_3)
               p3(1) = E_3
               p3(2) = px_3
               p3(3) = py_3
               p3(4) = pz_3
               v3 = p3/E_3
  
               x3 = x1+(v1-v3)*dt_coll
 
               n_jpsi_sur = n_jpsi_sur+1
               n_jpsi_pro = n_jpsi_pro+1

               if(sig_corr.eq.0)then
                  weight = 1.0/n_group/n_group
               elseif(gid(i_1).eq.gid(i_2))then
                  weight = 1.0/n_group
               else
                  weight = 1.0/n_group/n_group
               endif

               jpsi(1, n_jpsi_pro) = x3(2)
               jpsi(2, n_jpsi_pro) = x3(3)
               jpsi(3, n_jpsi_pro) = x3(4)
               jpsi(4, n_jpsi_pro) = -0.1
               jpsi(5, n_jpsi_pro) = p3(2)
               jpsi(6, n_jpsi_pro) = p3(3)
               jpsi(7, n_jpsi_pro) = p3(4)
               jpsi(8, n_jpsi_pro) = m_psi
               jpsi(9, n_jpsi_pro) = weight
 
               goto 3001
           
 3000          continue             

            enddo

 3001       continue

         enddo

      endsubroutine

      subroutine sr_jpsi_dissociation(jpsi, g, n_jpsi_pro, n_jpsi_sur,&
                                      n_g, m_psi, m_c, sigma, dt)

         implicit none

         !parameter variables=================================================
         ! array for the first species of particles
         ! The 8 components: 1-x, 2-y, 3-z, 4-t, 5-px, 6-py, 7-pz, 8-m
         real(8), intent(inout), dimension(9,200000) :: jpsi
         real(8), intent(inout), dimension(9,200000) :: g
  
         real(8), intent(in) :: sigma ! total cross section
         real(8), intent(in) :: dt    ! time step
         real(8), intent(in) :: m_psi, m_c
         integer(8), intent(inout) :: n_jpsi_pro, n_jpsi_sur
         integer(8), intent(in) :: n_g
  
         !====================================================================
  
         !temprary variables==================================================
         real(8) :: t  ! time 
  
         integer(8) :: i_1, i_2
         real(8) :: x_1, y_1, z_1, px_1, py_1, pz_1, m_1, E_1 ! particle_1
         real(8) :: x_2, y_2, z_2, px_2, py_2, pz_2, m_2, E_2 ! particle_2
  
         !_m for meeting (minimum distance in CMF)
  
         real(8) :: dis2, dis2_max, dis_min
         real(8) :: pi
         real(8) :: epsilon_v2
         real(8) :: s  ! square of the total energy in the center of mass frame
         real(8) :: ds ! ds= s-(m_1+m_2)^2, ds+4*m_1*m_2=s-(m_1-m_2)^2
  
         real(8) :: x1(4), x2(4)
         real(8) :: p1(4), p2(4), p(4)
         real(8) :: v1(4)
  
         real(8) :: dt_coll
         integer(8) :: sig_coll

         ! 0: update neither, 1: only update particle_1, 
         ! 2: only update particle_2, 3: update both
         !integer(8) :: sig_update
  
         !====================================================================
         !initialization======================================================
         pi = atan(1.0)*4.0
         dis2_max = sigma/pi + 4.0*dt*dt ! An upper boundary to collide 
         dis_min = 1.d-9 ! very small distance is not allowed
  
         epsilon_v2 = 1.d-14 ! to avoid v>=c
         t = 0.0
         !====================================================================
  
  
         !core code===========================================================
         do i_1 = 1, n_jpsi_pro
 
            m_1  = jpsi(8,i_1)
            if(m_1 .lt. 0.0) goto 3001
            x_1  = jpsi(1,i_1)
            y_1  = jpsi(2,i_1)
            z_1  = jpsi(3,i_1)
  
            px_1 = jpsi(5,i_1)
            py_1 = jpsi(6,i_1)
            pz_1 = jpsi(7,i_1)
  
            E_1  = sqrt(px_1**2+py_1**2+pz_1**2+m_1**2)
  
            p1(1)= E_1
            p1(2)= px_1
            p1(3)= py_1
            p1(4)= pz_1
            x1(1)= 0.0
            x1(2)= x_1
            x1(3)= y_1
            x1(4)= z_1
  
            v1 = p1/E_1
  
            do i_2 = 1, n_g
            
               x_2  = g(1,i_2)
               y_2  = g(2,i_2)
               z_2  = g(3,i_2)
  
               dis2 = (x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2)&
                     +(z_1-z_2)*(z_1-z_2)
               if(dis2 .gt. dis2_max) goto 3000
     
               px_2 = g(5,i_2)
               py_2 = g(6,i_2)
               pz_2 = g(7,i_2)
               m_2  = g(8,i_2)
  
               E_2 = sqrt(px_2**2+py_2**2+pz_2**2+m_2**2)
  
               p2(1)= E_2
               p2(2)= px_2
               p2(3)= py_2
               p2(4)= pz_2
               x2(1)= 0.0
               x2(2)= x_2
               x2(3)= y_2
               x2(4)= z_2
  
               call sr_coll_check (x1, x2, p1, p2, dt, sigma,&
                                   sig_coll, dt_coll)
               if(sig_coll .eq. 0) goto 3000
               p = p1+p2
               call sr_dot(p,p,s)
               ds = s-(m_c+m_c)**2 
               if( ds .le. 0.0) goto 3000

               ! Collision happens
 
               jpsi(8,i_2) = -m_psi
               n_jpsi_sur = n_jpsi_sur-1
 
               goto 3001
           
 3000          continue             

            enddo

 3001       continue

         enddo

      endsubroutine

      subroutine sr_free_streaming(particle, length, dt, x_up, x_low, &
            y_up, y_low, z_up, z_low, sig_cycl)
           implicit none
  
         !parameter variables=================================================
         ! array for particles
         ! The 8 components: 1-x, 2-y, 3-z, 4-t, 5-px, 6-py, 7-pz, 8-m
         real(8), intent(inout), dimension(9,200000) :: particle
  
         integer(8), intent(in) :: length ! number of particles in particle_1
         real(8), intent(in)    :: dt    ! time step
  
         ! for periodic boundary condition
         integer(8), intent(in) :: sig_cycl ! switch for periodic condition. 0-off, 1-on
         real(8), intent(in) :: x_up, y_up, z_up ! upper boundary
         real(8), intent(in) :: x_low, y_low, z_low ! lower boundary
         !====================================================================
  
         !temprary variables==================================================
         real(8) :: x, y, z
         real(8) :: px, py, pz, m, E
         real(8) :: vx, vy, vz
  
         integer(8) :: i
  
         real(8) :: length_x, length_y, length_z
         !====================================================================
  
         !streaming===========================================================
         do i=1, length
  
            ! initialization
            x  = particle(1,i)
            y  = particle(2,i)
            z  = particle(3,i)
            px = particle(5,i)
            py = particle(6,i)
            pz = particle(7,i)
            m  = particle(8,i)
  
            ! calculation velocity
            E  = sqrt(px*px+py*py+pz*pz+m*m)
            vx = px/E
            vy = py/E
            vz = pz/E
  
            ! calculating the new coordinates
            x  = x + vx*dt
            y  = y + vy*dt
            z  = z + vz*dt
  
            !checking the periodic boundary condition
            if(sig_cycl .eq. 1) then
               length_x = x_up - x_low
               length_y = y_up - y_low
               length_z = z_up - z_low
  
               if(x .gt. x_up)  x = x - length_x
               if(y .gt. y_up)  y = y - length_y
               if(z .gt. z_up)  z = z - length_z
  
               if(x .lt. x_low) x = x + length_x
               if(y .lt. y_low) y = y + length_y
               if(z .lt. z_low) z = z + length_z
  
            
               if(x .gt. x_up) print *, "Wrong x value ", x
            endif
  
            ! update the coordinates 
            particle(1,i) = x
            particle(2,i) = y
            particle(3,i) = z
  
         enddo
         !====================================================================
  
      endsubroutine sr_free_streaming




      !=================================================================
      ! initialize charm and anticharm
      !-----------------------------------
      subroutine sr_initialize_c(c, cbar, gid, n_pair, n_group, m_c,&
                                 dz, radius, sig_corr)
         implicit none
         real(8), intent(out) :: c(9, 200000)
         real(8), intent(out) :: cbar(9, 200000)
         real(8), intent(out) :: gid(200000)
         real(8), intent(in)  :: m_c, dz, radius
         integer(8), intent(in) :: n_pair, n_group, sig_corr

         integer(8) :: i, n_c
         real(8) :: x, y, z, px, py, pz, E_c
         real(8) :: r, phi
         real(8) :: pi
         real(8) :: rand2, f, pt
         pi = atan(1.0)*4.0

         n_c = n_pair*n_group

         do i=1, n_c

            gid(i) = (i-1)/n_group+1

            phi = 2.0*pi*rand()
            r = radius*(rand())**(1/2.0) 
  
            x = r*cos(phi)
            y = r*sin(phi)
            z = dz*(rand()-0.5) 
      
          
 12         pt = 6.0*rand()
            f = pt*(1+(pt/6)**2)/(1+pt/3.7)**12/&
                (1+exp(0.9-2.0*pt))
            rand2 = 0.06*rand() 
            if (rand2 .gt. f) goto 12  

            phi = 2.0*pi*rand() 

            px = pt*cos(phi)
            py = pt*sin(phi)
            pz = sqrt(m_c**2+pt**2)*sinh(rand()-0.5) 
      


            E_c = sqrt(m_c**2+pt**2)
            c(1,i) = x 
            c(2,i) = y
            c(3,i) = z 
            c(4,i) = 1.0/(2.0*E_c) 
            c(5,i) = px
            c(6,i) = py 
            c(7,i) = pz 
            c(8,i) = m_c 

            if (sig_corr .eq. 1) then
               cbar(1,i) =  c(1,i)
               cbar(2,i) =  c(2,i)
               cbar(3,i) =  c(3,i)
               cbar(4,i) =  c(4,i)
               cbar(5,i) = -c(5,i)
               cbar(6,i) = -c(6,i)
               cbar(7,i) =  c(7,i)
               cbar(8,i) =  c(8,i)
            else
               phi = 2.0*pi*rand()
               r = radius*(rand())**(1/2.0) 
     
               x = r*cos(phi)
               y = r*sin(phi)
               z = dz*(rand()-0.5) 
         
             
 13            pt= 6.0*rand()
               f = pt*(1+(pt/6)**2)/(1+pt/3.7)**12/&
                   (1+exp(0.9-2.0*pt))
               rand2 = 0.06*rand() 
               if (rand2 .gt. f) goto 13  
   
               phi = 2.0*pi*rand() 
   
               px = pt*cos(phi)
               py = pt*sin(phi)
               pz = sqrt(m_c**2+pt**2)*sinh(rand()-0.5) 
         
               E_c = sqrt(m_c**2+pt**2)
               cbar(1,i) = x 
               cbar(2,i) = y
               cbar(3,i) = z 
               cbar(4,i) = 1.0/(2.0*E_c) 
               cbar(5,i) = px
               cbar(6,i) = py 
               cbar(7,i) = pz 
               cbar(8,i) = m_c 
            endif
         end do
      end subroutine

      subroutine sr_jpsi_sur (jpsi, n_jpsi_pro, res)
         real(8), intent(in) :: jpsi(9,200000)
         integer(8), intent(in) :: n_jpsi_pro
         real(8), intent(out):: res

         integer(8) :: i
         res = 0.0
         do i=1, n_jpsi_pro
            if(jpsi(8,i).gt.0.0) then
               res = res+jpsi(9,i)
            endif
         enddo
      endsubroutine

      subroutine sr_jpsi_pro (jpsi, n_jpsi_pro, res)
         real(8), intent(in) :: jpsi(9,200000)
         integer(8), intent(in) :: n_jpsi_pro
         real(8), intent(out):: res

         integer(8) :: i
         res = 0.0
         do i=1, n_jpsi_pro
            res = res+jpsi(9,i)
         enddo
      endsubroutine
           
