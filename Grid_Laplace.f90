program Grid_Laplace

    implicit none
   
!*************** Definition of Variables ***************!   
        integer, parameter :: IM = 6, JM = 6, kmax = IM*JM, tot_it = 2000

        real(kind = 8), parameter :: Lx = 2, Ly = 1, dx = Lx/(IM-1), dy = Ly/(JM-1), &
                                     tol = 1e-6, w = .8

        integer :: i, j, k, iter_x, last_iter_x, iter_y, last_iter_y

        real(kind = 8) :: pi, phi, & 
                          x_ksi, x_eta, x_ksi_ksi, x_eta_eta, x_ksi_eta, &
                          y_ksi, y_eta, y_ksi_ksi, y_eta_eta, y_ksi_eta, &
                          g_1_Cova_x, g_1_Cova_y, g_2_Cova_x, g_2_Cova_y, &
                          g_11_Cova, g_22_Cova, g_12_Cova, &
                          g_1_Contra_x, g_1_Contra_y, g_2_Contra_x, g_2_Contra_y, &
                          start, finish, &
                          norm_x, rms_x, min_res_x, &
                          norm_y, rms_y, min_res_y, &
                          wall_clock_time

        real(kind = 8), dimension(kmax) :: x, x_old, y, y_old, Jac, &
                                           g_11_Contra, g_22_Contra, g_12_Contra, &
                                           d1, d2, d3, d4, d5, &
                                           residual_x, residual_y

        real(kind=8), dimension(tot_it) :: norm_array_x, norm_array_y
                                        
!******************* Initialization ********************!
                do i = 1, IM
                    do j = 1, JM
                        k = ((i-1)*JM) + j

                        Jac(k) = 0.d0

                        d1(k) = 0.d0
                        d2(k) = 1.d0
                        d3(k) = 0.d0
                        d4(k) = 0.d0
                        d5(k) = 0.d0

                        x(k) = 0.d0
                        x_old(k) = 0.d0
                        y(k) = 0.d0
                        y_old(k) = 0.d0
                        residual_x(k) = 0.d0
                        residual_y(k) = 0.d0
                    enddo
                enddo

            iter_x = 0					
            norm_x = 1					
            rms_x = 1

            iter_y = 0					
            norm_y = 1					
            rms_y = 1

!******************* Grid Generation *******************! Generation of Boundary Nodes
            phi = 0.d0
            pi = 4.d0*atan(1.d0) 
            phi=phi*pi/180.d0
                do i = 1, IM
                    do j = 1, JM
                        k = ((i-1)*JM) + j
                            if (i.eq.1.or.i.eq.IM.or.j.eq.1.or.j.eq.JM) then 
                                x(k)=(i*dx*dcos(phi))-(j*dy*dsin(phi))-dx
                                y(k)=(i*dx*dsin(phi))+(j*dy*dcos(phi))-dy
                            endif
                    enddo
                enddo

!**************** Bilinear Interpolation ***************! Generation of Internal Nodes
            do i = 2, (IM-1)
                do j = 2, (JM-1)
                    k = ((i-1)*JM) + j
                        x(k)= (((real(i)-real(IM))/(1.d0-real(IM)))*((real(j)-real(JM))/(1.d0-real(JM)))*x(1)) + &
                              (((real(i)-1.d0)/(real(IM)-1.d0))*((real(j)-real(JM))/(1.d0-real(JM)))*x(kmax-JM+1)) + &
                              (((real(i)-real(IM))/(1.d0-real(IM)))*((real(j)-1.d0)/(real(JM)-1.d0))*x(JM)) + &
                              (((real(i)-1.d0)/(real(IM)-1.d0))*((real(j)-1.d0)/(real(JM)-1.d0))*x(kmax))

                        y(k)= (((real(i)-real(IM))/(1.d0-real(IM)))*((real(j)-real(JM))/(1.d0-real(JM)))*y(1)) + &
                              (((real(i)-1.d0)/(real(IM)-1.d0))*((real(j)-real(JM))/(1.d0-real(JM)))*y(kmax-JM+1)) + &
                              (((real(i)-real(IM))/(1.d0-real(IM)))*((real(j)-1.d0)/(real(JM)-1.d0))*y(JM)) + &
                              (((real(i)-1.d0)/(real(IM)-1.d0))*((real(j)-1.d0)/(real(JM)-1.d0))*y(kmax))
                enddo
            enddo

            x(9) = 4
            y(9) = 4

!**************** Initial Grid Printout ****************!
            open(unit = 1, file = "Initial_Grid.txt")
            
                do i = 1, IM
                    do j = 1, JM
                        k = ((i-1)*JM) + j
                            write(1,*) x(k), y(k)
                    enddo
                        write(1,*)
                enddo

                do j = 1, JM    ! The reverse loop, is for the grid visualization
                    do i = 1, IM
                        k = ((i-1)*JM) + j
                            write(1,*) x(k), y(k)
                    enddo
                        write(1,*)
                enddo

            close(unit=1)

!*************** Boundary Nodes Printout ***************!
            open(unit = 2, file = "Boundary_Nodes.txt")

                i = 1   ! Left Boundary
                    do j = 1, JM
                        k = ((i-1)*JM) + j
                            write(2,*) x(k), y(k)
                    enddo

                j = JM  ! Top Boundary
                    do i = 2, IM
                        k = ((i-1)*JM) + j
                            write(2,*) x(k), y(k)
                    enddo

                i = IM  !Right Boundary
                    do j = (JM-1), 1, -1
                        k = ((i-1)*JM)+j
                            write(2,*) x(k), y(k)
                    enddo

                j = 1   !Bottom Boundary
                    do i = (IM-1), 1, -1
                        k = ((i-1)*JM)+j
                            write(2,*) x(k), y(k)
                    enddo
            
            close(unit = 2)

!**************** Derivatives & Metrices ***************!
            open(unit = 3, file = "Jacobian.txt")
            open(unit = 4, file = "Contravariant_Metrices.txt")

                do i = 2, (IM-1)
                    do j = 2, (JM-1)
                        k = ((i-1)*JM) + j

!************** Derivatives (x-Direction) **************!
                            x_ksi = (x(k+JM) - x(k-JM))/2.d0
                            x_eta = (x(k+1) - x(k-1))/2.d0
                            x_ksi_ksi = x(k+JM) - (2.d0*x(k)) + x(k-JM)
                            x_eta_eta = x(k+1) - (2.d0*x(k)) + x(k-1)
                            x_ksi_eta = (1.d0/4.d0)*(x(k+JM+1) - x(k-JM+1) - x(k+JM-1) + x(k-JM-1))

!************** Derivatives (y-Direction) **************!
                            y_ksi = (y(k+JM) - y(k-JM))/2.d0
                            y_eta = (y(k+1) - y(k-1))/2.d0
                            y_ksi_ksi = y(k+JM) - (2.d0*y(k)) + y(k-JM)
                            y_eta_eta = y(k+1) - (2.d0*y(k)) + y(k-1)
                            y_ksi_eta = (1.d0/4.d0)*(y(k+JM+1) - y(k-JM+1) - y(k+JM-1) + y(k-JM-1))

!*********************** Jacobian **********************!
                            Jac(k) = (x_ksi*y_eta) - (x_eta*y_ksi)
                                write(3,*) k, Jac(k)

!************* 1st Order Covariant Metrices ************!
!************************** g1 *************************!
                            g_1_Cova_x = x_ksi
                            g_1_Cova_y = y_ksi

!************************** g2 *************************!
                            g_2_Cova_x = x_eta
                            g_2_Cova_y = y_eta

!************* 2nd Order Covariant Metrices ************!
!************************** g11 ************************!
                            g_11_Cova = (g_1_Cova_x**2) + (g_1_Cova_y**2)

!************************** g22 ************************!
                            g_22_Cova = (g_2_Cova_x**2) + (g_2_Cova_y**2)

!************************** g12 ************************!
                            g_12_Cova = (g_1_Cova_x*g_2_Cova_x) + (g_1_Cova_y*g_2_Cova_y)

!*********** 1st Order Contravariant Metrices **********!
!************************** g1 *************************!
                            g_1_Contra_x = y_eta/Jac(k)
                            g_1_Contra_y = - (x_eta/Jac(k))

!************************** g2 *************************!
                            g_2_Contra_x = - (y_ksi/Jac(k))
                            g_2_Contra_y = x_ksi/Jac(k)

!*********** 2nd Order Contravariant Metrices **********!
!************************** g11 ************************!
                            g_11_Contra(k) = (g_1_Contra_x**2) + (g_1_Contra_y**2)

!************************** g22 ************************!
                            g_22_Contra(k) = (g_2_Contra_x**2) + (g_2_Contra_y**2)

!************************** g12 ************************!
                            g_12_Contra(k) = (g_1_Contra_x*g_2_Contra_x) + (g_1_Contra_y*g_2_Contra_y)
                                write(4,*) k, &
                                           g_1_Contra_x, g_1_Contra_y, &
                                           g_2_Contra_x, g_2_Contra_y, &
                                           g_11_Contra(k), g_22_Contra(k), g_12_Contra(k)

!********************** Diagonals **********************!
                            d1(k) = ((1.d0/2.d0)*g_12_Contra(k)) - g_22_Contra(k) 
                            d2(k) = 2.d0*(g_11_Contra(k) + g_22_Contra(k))
                            d3(k) = ((1.d0/2.d0)*g_12_Contra(k)) - g_22_Contra(k)  
                            d4(k) = - g_11_Contra(k) - ((1.d0/2.d0)*g_12_Contra(k))
                            d5(k) = - g_11_Contra(k) - ((1.d0/2.d0)*g_12_Contra(k)) 

                    enddo
                        write(3,*)
                        write(4,*)
                enddo

            close(unit = 3)
            close(unit = 4)

!************************* G-S *************************!
            call cpu_time(start)

!****************** G-S (x-Coordinate) *****************!
                do while ((rms_x > tol).and.(iter_x < tot_it))
                    iter_x = iter_x + 1
                    x_old = x

                        do i = 2, (IM-1)
                            do j = 2, (JM-1)
                                k = ((i-1)*JM) + j

                                x(k) = ((1-w)*x_old(k)) + ((w/d2(k))*(- (d1(k)*x(k-1)) &
                                                                      - (d3(k)*x_old(k+1)) &
                                                                      - (d4(k)*x(k-JM)) &
                                                                      - (d5(k)*x_old(k+JM))))

                            enddo
                        enddo

!*********************** Residual **********************!
                            do i = 2, (IM-1)
                                do j = 2, (JM-1)
                                    k = ((i-1)*JM) + j

                                    residual_x(k) = - &
                                                    ((d1(k)*x(k-1)) + &
                                                    (d2(k)*x(k)) + &
                                                    (d3(k)*x(k+1)) + &
                                                    (d4(k)*x(k-JM)) + &
                                                    (d5(k)*x(k+JM)))

                                enddo
                            enddo

                                norm_x = norm2(residual_x)
                                rms_x = norm_x/(sqrt(real(kmax)))
                                norm_array_x(iter_x) = norm_x

                enddo

                    last_iter_x = iter_x
                    min_res_x = minval(norm_array_x)

!****************** G-S (y-Coordinate) *****************!
                do while ((rms_y > tol).and.(iter_y < tot_it))
                    iter_y = iter_y + 1
                    y_old = y
    
                        do i = 2, (IM-1)
                            do j = 2, (JM-1)
                                k = ((i-1)*JM) + j
    
                                y(k) = ((1-w)*y_old(k)) + ((w/d2(k))*(- (d1(k)*y(k-1)) &
                                                                      - (d3(k)*y_old(k+1)) &
                                                                      - (d4(k)*y(k-JM)) &
                                                                      - (d5(k)*y_old(k+JM))))
    
                            enddo
                        enddo
    
!*********************** Residual **********************!
                            do i = 2, (IM-1)
                                do j = 2, (JM-1)
                                    k = ((i-1)*JM) + j
    
                                    residual_y(k) = - &
                                                    ((d1(k)*y(k-1)) + &
                                                    (d2(k)*y(k)) + &
                                                    (d3(k)*y(k+1)) + &
                                                    (d4(k)*y(k-JM)) + &
                                                    (d5(k)*y(k+JM)))
    
                                enddo
                            enddo
    
                                norm_y = norm2(residual_y)
                                rms_y = norm_y/(sqrt(real(kmax)))
                                norm_array_y(iter_y) = norm_y
    
                enddo
    
                    last_iter_y = iter_y
                    min_res_y = minval(norm_array_y)

                    wall_clock_time = start - finish

            call cpu_time(finish)
  
!*********************** Printout **********************!
            open(unit = 5, file = "Final_Grid.txt")          
            
                do i = 1, IM
                    do j = 1, JM
                        k = ((i-1)*JM) + j
                            write(5,*) x(k), y(k)
                    enddo
                        write(5,*)
                enddo

                do j = 1, JM    ! The reverse loop, is for the grid visualization
                    do i = 1, IM
                        k = ((i-1)*JM) + j
                            write(5,*) x(k), y(k)
                    enddo
                        write(5,*)
                enddo

            close(unit = 5)

            open(unit = 6, file = "Nodes_min_max.txt")

                write(6,*) "Min_x:",minval(x), "|", "Max_x:", maxval(x)
                write(6,*) "Min_y:",minval(y), "|", "Max_y:", maxval(y)

            close(unit = 6)

            open(unit = 7, file = "Norm_x.txt")

                do iter_x = 1, last_iter_x
                    write(7,*) iter_x, norm_array_x(iter_x)
                enddo

            close(unit = 7)

            open(unit = 8, file = "Norm_y.txt")

                do iter_y = 1, last_iter_y
                    write(8,*) iter_y, norm_array_y(iter_y)
                enddo

            close(unit = 8)

            print *, "The CPU needed-time of the algorithm was:", wall_clock_time, "(sec)"
            print *, ("")
            print *, "The total number of iterations for the abscissa &
                     & (x-axis value) was:", last_iter_x
            print *, ("")
            print *, "The total number of iterations for the ordinate & 
                     & (y-axis value) was:", last_iter_y
            print *, ("")
            print *, "The minimum residual of the final iteration for &
                     & the abscissa (x-axis value) was:", min_res_x
            print *, ("")
            print *, "The minimum residual of the final iteration for & 
                     & the ordinate (y-axis value) was:", min_res_y
            print *, ("")
            print *, "The requested results have been written on the & 
                     & corresponding files. Thank you!"
            print *, ("")

end program Grid_Laplace