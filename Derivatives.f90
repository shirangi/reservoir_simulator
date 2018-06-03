!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! 3D oil reservoir Simulator project 
!! by Mehrdad Gharib Shirangi From Spring 2010 until Fall 2010
!! mehrdad.ghsh@gmail.com
!! Copy Right by Mehrdad Gharib Shirangi
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Module Derivative

! contains transmissibility functions
use globalvar
use Trans

implicit none

!real :: Press

contains
   real*8 function DTrxa_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 ::  delx_avg , Press
   real*8 :: denomenator, press_term, Dip_cor
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + 1
      delx_avg = (Deltax(i+1) + Deltax(i)) / 2.0
      denomenator = 1 / 2.0 * ( Deltax(i+1) / ( kx(ind+1) * Deltay(j) * Deltaz(ind+1) ) + Deltax(i)/( kx(ind)*Deltay(j)*Deltaz(ind) ) )         
      Dip_cor = delx_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + delx_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then !( Potential(ind,m) > Potential(N_ind,m) ) then ! upstream block
         if (wrt=='a') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) + mu_der(ind ,'w') / mu_w(ind) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) + mu_der(N_ind ,'w') / mu_w(N_ind) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTrxa_DP = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrxa_DP
!******************************************************************************************* 
   real*8 function DTrxb_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 ::  delx_avg , Press
   real*8 :: denomenator, press_term, Dip_cor
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'b' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - 1
      delx_avg = (Deltax(i-1) + Deltax(i))/2.0
      denomenator = 1 / 2.0 * ( Deltax(i-1) / ( kx(ind-1) * Deltay(j) * Deltaz(ind-1) ) + Deltax(i)/( kx(ind)*Deltay(j)*Deltaz(ind) ) )         
      Dip_cor = delx_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + delx_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTrxb_DP = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrxb_DP   
!*******************************************************************************************
   real*8 function DTryb_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: dely_avg , Press
   real*8 :: denomenator, press_term, Dip_cor   
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'b' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - Nx
      dely_avg = (Deltay(j-1) + Deltay(j))/2.0
      denomenator = 1 / 2.0 * ( Deltay(j-1) / ( ky(N_ind) * Deltax(i) * Deltaz(N_ind) ) + Deltay(j)/( ky(ind) * Deltax(i) * Deltaz(ind) ) )         
      Dip_cor = dely_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + dely_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTryb_DP = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTryb_DP   
!*******************************************************************************************
   real*8 function DTrya_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: dely_avg, Press
   real*8 :: denomenator, press_term, Dip_cor      
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + Nx
      dely_avg = (Deltay(j+1) + Deltay(j))/2.0
      denomenator = 1 / 2.0 * ( Deltay(j+1) / ( ky(N_ind) * Deltax(i) * Deltaz(N_ind) ) + Deltay(j)/( ky(ind) * Deltax(i) * Deltaz(ind) ) )         
      Dip_cor = dely_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + dely_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0.0 ) then ! upstream block
         if (wrt=='a') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTrya_DP = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrya_DP
   
   
!*******************************************************************************************

   real*8 function DTrzb_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: denomenator, press_term , Press
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - Nx * Ny
      denomenator = 1 / 2.0 * ( Deltaz(N_ind) / ( kz(N_ind) * Deltax(i) * Deltay(j) ) + Deltaz(ind)/( kz(ind) * Deltax(i) * Deltay(j) ) )         
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTrzb_DP = betac * press_term /denomenator      
      
   Endfunction DTrzb_DP
   
!*******************************************************************************************

   real*8 function DTrza_DP(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 ::delz_avg , Press
   real*8 :: denomenator, press_term    
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + Nx * Ny
      denomenator = 1 / 2.0 * ( Deltaz(N_ind) / ( kz(N_ind) * Deltax(i) * Deltay(j) ) + Deltaz(ind)/( kz(ind) * Deltax(i) * Deltay(j) ) )         
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='a') then
            press_term = 0.0
         Else
            if ( m == 'w') then
			    press_term = - k_rw(ind) / ( mu_w (ind) * B_w(ind ) ) * ( B_der(ind ,'w') / B_w(ind ) )               
            Else
				press_term = - k_ro(ind) / ( mu_o (ind) * B_o(ind ) ) * ( B_der(ind ,'o') / B_o(ind ) + mu_der(ind , 'o') / mu_o(ind) )               
            EndIf                    
         EndIf
      else
         if (wrt=='c') then
            press_term = 0.0
         Else
             if ( m == 'w') then
			    press_term = - k_rw( N_ind ) / ( mu_w (N_ind) * B_w(N_ind ) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) )                               
             Else
				press_term = - k_ro( N_ind) / ( mu_o ( N_ind) * B_o( N_ind ) ) * ( B_der( N_ind ,'o') / B_o( N_ind ) + mu_der( N_ind , 'o') / mu_o( N_ind) )                               
             EndIf                             
         Endif 
      EndIf
      
      DTrza_DP = betac * press_term /denomenator      
      
   Endfunction DTrza_DP
   
!*******************************************************************************************
     real*8 function DTrxa_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: delx_avg , Press
   real*8 :: denomenator, press_term, Dip_cor         
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + 1
      delx_avg = (Deltax(i+1) + Deltax(i))/2.0
      denomenator = 1 / 2.0 * ( Deltax(i+1) / ( kx(ind+1) * Deltay(j) * Deltaz(ind+1) ) + Deltax(i)/( kx(ind)*Deltay(j)*Deltaz(ind) ) )         
      Dip_cor = delx_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + delx_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='a') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) ) + k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) + mu_der(ind ,'w') / mu_w(ind) ) * PCO_der(ind)
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) + mu_der(N_ind ,'w') / mu_w(N_ind) ) * PCO_der(N_ind)                            
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTrxa_DS = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrxa_DS
!*******************************************************************************************   
       real*8 function DTrxb_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: delx_avg , Press
   real*8 :: denomenator, press_term, Dip_cor         
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'b' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - 1
      delx_avg = (Deltax(i-1) + Deltax(i))/2.0
      denomenator = 1 / 2.0 * ( Deltax(i-1) / ( kx(N_ind) * Deltay(j) * Deltaz(N_ind) ) + Deltax(i)/( kx(ind)*Deltay(j)*Deltaz(ind) ) )         
      Dip_cor = delx_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + delx_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) ) + k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) ) * PCO_der(ind)                
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) ) * PCO_der(N_ind)                                                        
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTrxb_DS = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrxb_DS
!*******************************************************************************************  
  
   real*8 function DTryb_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 ::dely_avg , Press
   real*8 :: denomenator, press_term, Dip_cor         
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'b' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - Nx
      dely_avg = (Deltay(j-1) + Deltay(j))/2.0
      denomenator = 1 / 2.0 * ( Deltay(j-1) / ( ky(N_ind) * Deltax(i) * Deltaz(N_ind) ) + Deltay(j)/( ky(ind) * Deltax(i) * Deltaz(ind) ) )         
      Dip_cor = dely_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + dely_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) )+ k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) ) * PCO_der(ind)                 
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) ) * PCO_der(N_ind)                                                        
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTryb_DS = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTryb_DS

!*******************************************************************************************

   real*8 function DTrya_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: dely_avg , Press
   real*8 :: denomenator, press_term, Dip_cor         
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + Nx
      dely_avg = (Deltay(j+1) + Deltay(j) )/2.0
      denomenator = 1 / 2.0 * ( Deltay(j+1) / ( ky(N_ind) * Deltax(i) * Deltaz(N_ind) ) + Deltay(j)/( ky(ind) * Deltax(i) * Deltaz(ind) ) )         
      Dip_cor = dely_avg ** 2 / ((Zblock(ind)-Zblock(N_ind)) ** 2 + dely_avg ** 2)
      
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='a') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) ) + k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) ) * PCO_der(ind)                
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) ) * PCO_der(N_ind)                                                        
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTrya_DS = betac * Dip_cor * press_term /denomenator      
      
   Endfunction DTrya_DS
!*******************************************************************************************
  
   real*8 function DTrzb_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: denomenator, press_term , Press
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'b' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind - Nx * Ny
      denomenator = 1 / 2.0 * ( Deltaz(N_ind) / ( kz(N_ind) * Deltax(i) * Deltay(j) ) + Deltaz(ind)/( kz(ind) * Deltax(i) * Deltay(j) ) )         
            
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='b') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) ) + k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) ) * PCO_der(ind)                 
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) ) * PCO_der(N_ind)                                                        
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTrzb_DS = betac * press_term /denomenator      
      
   Endfunction DTrzb_DS
!*******************************************************************************************

  
   real*8 function DTrza_DS(i,j,k,m,wrt) ! the value of the derivative of x transmissibility between a gridblock and a next grid 
   integer :: ind , N_ind ! Neighbor grid index
   integer :: i,j,k
   real*8 :: denomenator, press_term , Press
   Character*1 :: m ! 'o' or 'w' : oil or water
   Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
      
      ind = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_ind = ind + Nx * Ny
      denomenator = 1 / 2.0 * ( Deltaz(N_ind) / ( kz(N_ind) * Deltax(i) * Deltay(j) ) + Deltaz(ind)/( kz(ind) * Deltax(i) * Deltay(j) ) )         
            
      if ( D_Poten( ind , N_ind , m ) > 0 ) then ! upstream block
         if (wrt=='a') then
            press_term = 0            
         Else
             if ( m == 'w') then
				press_term = DkrwdS( ind ) / ( mu_w ( ind ) * B_w( ind ) ) + k_rw(ind) / ( mu_w ( ind) * B_w(ind) ) * ( B_der(ind ,'w') / B_w(ind ) ) * PCO_der(ind)                
             Else
				press_term = DkrodS( ind ) / ( mu_o ( ind ) * B_o( ind ) )                                 
             EndIf                             
         EndIf
      else
         if (wrt=='c') then
            press_term = 0
         Else
            if ( m == 'w') then
				press_term = DkrwdS( N_ind ) / ( mu_w ( N_ind ) * B_w( N_ind ) ) + k_rw(N_ind) / ( mu_w ( N_ind) * B_w(N_ind) ) * ( B_der(N_ind ,'w') / B_w(N_ind ) ) * PCO_der(N_ind)                                                        
            Else
				press_term = DkrodS( N_ind ) / ( mu_o ( N_ind ) * B_o( N_ind ) )                                             
            EndIf                                               
         Endif 
      EndIf
      
      DTrza_DS = betac * press_term /denomenator      
      
   Endfunction DTrza_DS 
!*******************************************************************************************
!*******************************************************************************************

  real*8 function B_der(ind , m )
    integer :: ind , n2 , m2 , i4   ! Neighbor grid index
    Character*1 :: m ! 'o' or 'w' : oil or water
    real*8 :: p2 , X
    
           
           if (m=='o') then
			  p2 = pressure(ind)           
              n2 = N_Bo
              m2 = 1
              do while (abs(n2-m2) .NE. 1)
                i4 = FLOOR(m2 + (n2-m2+1)/2.0)
                if (p2 < Bo_p(i4,1) )  then
                    n2 = FLOOR(m2+(n2-m2+1)/2.0);
                else
                    m2 = FLOOR(m2+(n2-m2+1)/2.0);
                end if
              end do
              B_der = (Bo_p(n2,2)-Bo_p(m2,2))/(Bo_p(n2,1)-Bo_p(m2,1))
!			   B_der = -5.0d-5
           Else
              p2 = W_press(ind)
              !B_der = - PVTW(2) / ( ( 1 + PVTW(3) * (p2 - PVTW(1) ) )**2 ) * PVTW(3)
              X =  PVTW(3) * (p2 - PVTW(1) )
              B_der = - PVTW(2) / ( ( 1 + X + ( X ** 2 )/2.0 )**2 ) * ( PVTW(3) + ( PVTW(3) ** 2 ) * (p2 - PVTW(1) ) )
              
           EndIf
  Endfunction B_der
  
  !*******************************************************************************************
  
  real*8 function mu_der(ind , m )
    integer :: ind ! Neighbor grid index
    Character*1 :: m ! 'o' or 'w' : oil or water
    real*8 :: p1
    Integer :: m1,n1 , i5
    
    p1 = pressure(ind)   
    if (m=='w') then
        mu_der = 0.0
    else
        n1 = N_mup
        m1 = 1
        do while (abs(n1-m1) .NE. 1)
            i5 = floor(m1 + (n1-m1+1)/2.0)
            if (p1 < mu_p(i5,1) )  then
                n1 = floor(m1+(n1-m1+1)/2.0);
            else
                m1 = floor(m1+(n1-m1+1)/2.0);
            end if
        end do
        mu_der = (mu_p(n1,2)-mu_p(m1,2))/(mu_p(n1,1)-mu_p(m1,1))
    Endif

  Endfunction mu_der
    
  !*******************************************************************************************
  
  real*8 function gama_der(ind , m)
    integer :: ind ! Ind is the block number which we want to have the derivative wrt the pressure of that block
    Character*1 :: m ! 'o' or 'w' : oil or water
    !Character*1 :: wrt ! 'a' , 'c'  ! wrt to after , or before or itself
     
    If ( m == 'o' ) then 
       gama_der = - den_o / ( 144.0 * 2.0 * ( B_o(ind) ) **2 ) *  B_der(ind , 'o' ) 
    Else !( m == 'w' ) 
       gama_der = - den_w / ( 144.0 * 2.0 * ( B_w(ind) ) **2 ) *  B_der(ind , 'w' ) 
    EndIf
    
  Endfunction gama_der
  !******************************************************************************************* 

  real*8 function Dkr_dS(index , m1)   ! this function calculates the Derivative of Kr wrt the saturation of the block

    Implicit None
    !
    integer :: index ! Ind is the block number which we want to have the derivative wrt the saturation of that block
    Character*1 :: m1 ! 'o' or 'w' : oil or water 
    real*8 :: swater
    Integer :: m ,n 

        swater = sw(index)
        n = n_Kr
        m = 1
        do while (abs(n-m) .NE. 1)
            if ( swater < Sw_mat(floor(m+(n-m+1)/2.0)))  then
                n = floor(m+(n-m+1)/2.0);
            else
                m = floor(m+(n-m+1)/2.0);
            end if
        end do
        if  ( m1 == 'o' ) then
			if ((swater > Sw_mat(n_Kr - 1) ) .OR. (swater .LT. Sw_mat(1) ) ) then
				Dkr_dS = 0.0
			Else        
				Dkr_dS = (kro(n)-kro(m))/(Sw_mat(n)-Sw_mat(m))
			EndIf            
        else
			if ( (swater > 1.0 ) .OR. (swater .LT. Sw_mat(1) ) ) then
				Dkr_dS = 0.0
			Else        
				Dkr_dS = (krw(n)-krw(m))/(Sw_mat(n)-Sw_mat(m))
			EndIf            
        EndIf     

  Endfunction Dkr_dS
  
  !******************************************************************************************* 
  real*8 function PCO_der(ind)
    integer :: ind 
    real*8 :: Swater
    Integer :: m,n
        
        Swater = sw(ind)

		if ( Swater .LE. Sw_mat(1) ) then
			n = 2
			m = 1
		elseif ( Swater .GE. 1.0   ) then 
			n = n_Kr
			m = n_Kr - 1
		else		
			n = n_Kr
			m=1
			do while (abs(n-m) .NE. 1)
				if (Swater < Sw_mat( floor(m+(n-m+1)/2.0) ) ) then
					n = floor(m+(n-m+1)/2.0);
				else
					m = floor(m+(n-m+1)/2.0);
				end if
			end do		
		EndIf
        PCO_der = (Pcow(n) - Pcow(m))/(Sw_mat(n) - Sw_mat(m))
	   ! PCO_der = - 3 * 34.98 * Swater ** 2 + 2 * 83.96 * Swater - 67.17
        

  Endfunction PCO_der
  !*******************************************************************************************
!*******************************************************************************************

  real*8 function B_der_P(p , m )
    integer :: ind , n2 , m2 , i4   ! Neighbor grid index
    Character*1 :: m ! 'o' or 'w' : oil or water
    real*8 :: p , p2 , X
    
           p2 = p           
           if (m=='o') then
              n2 = N_Bo
              m2 = 1
              do while (abs(n2-m2) .NE. 1)
                i4 = FLOOR(m2 + (n2-m2+1)/2.0)
                if (p2 < Bo_p(i4,1) )  then
                    n2 = FLOOR(m2+(n2-m2+1)/2.0);
                else
                    m2 = FLOOR(m2+(n2-m2+1)/2.0);
                end if
              end do
              B_der_P = (Bo_p(n2,2)-Bo_p(m2,2))/(Bo_p(n2,1)-Bo_p(m2,1))
!			   B_der = -5.0d-5
           Else
              !B_der = - PVTW(2) / ( ( 1 + PVTW(3) * (p2 - PVTW(1) ) )**2 ) * PVTW(3)
              X =  PVTW(3) * (p2 - PVTW(1) )
              B_der_P = - PVTW(2) / ( ( 1 + X + ( X ** 2 )/2.0 )**2 ) * ( PVTW(3) + ( PVTW(3) ** 2 ) * (p2 - PVTW(1) ) )              
           EndIf
           
  Endfunction B_der_P
  
  !*******************************************************************************************
  !******************************************************************************************* 
  real*8 function PCO_der_Sw(Swater)
    integer :: ind 
    real*8 :: Swater
    Integer :: m,n        

		if ( Swater .LE. Sw_mat(1) ) then
			n = 2
			m = 1
		elseif ( Swater .GE. 1.0   ) then 
			n = n_Kr
			m = n_Kr - 1
		else		
			n = n_Kr
			m=1
			do while (abs(n-m) .NE. 1)
				if (Swater < Sw_mat( floor(m+(n-m+1)/2.0) ) ) then
					n = floor(m+(n-m+1)/2.0);
				else
					m = floor(m+(n-m+1)/2.0);
				end if
			end do		
		EndIf
        PCO_der_Sw = (Pcow(n) - Pcow(m))/(Sw_mat(n) - Sw_mat(m))
	   ! PCO_der = - 3 * 34.98 * Swater ** 2 + 2 * 83.96 * Swater - 67.17
        

  Endfunction PCO_der_Sw
  !*******************************************************************************************
  
  
  EndModule Derivative