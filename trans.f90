!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! 3D oil reservoir Simulator project 
!! by Mehrdad Gharib Shirangi From Spring 2010 until Fall 2010
!! mehrdad.ghsh@gmail.com
!! Copy Right by Mehrdad Gharib Shirangi
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Module Trans

! contains transmissibility functions
use globalvar
!!use  Derivative
!************************************************************************************************************************************************************   
!                                                     Explanations
 

!************************************************************************************************************************************************************   
implicit none

real*8 :: Dip_cor, & ! dip correction which is cos(theta)
          delax_avg , delay_avg ,  denomenator , press_term ! pressure dependent term
!!integer:: i1 , j1 



contains
   real*8 function Transxb(i,j,k,m) ! the value of the transmissibility between a gridblock and a grid before that ! gridblock (i-1,j,k)\\
   integer :: Index , N_Index
   integer :: i,j,k
   Character*1 :: m
   real*8 :: Dip_cor, delax_avg, denomenator, press_term 
   real*8 :: Press
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i  
      N_Index = Index - 1
      delax_avg = ( Deltax(i-1) + Deltax(i)) / 2.0
      Dip_cor = delax_avg ** 2 / ((Zblock(Index)-Zblock(Index-1)) ** 2 + delax_avg ** 2)
      denomenator = 1 / 2.0 * ( Deltax(i-1) / ( kx(Index-1) * Deltay(j) * Deltaz(Index-1) ) + Deltax(i)/( kx(Index)*Deltay(j)*Deltaz(Index) ) )
      
      if  ( D_Poten(Index , Index - 1 , m ) > 0 ) then ! (Potential(Index,m) > Potential(Index - 1,m) ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transxb = betac * Dip_cor * press_term /denomenator
   endfunction Transxb
!***************************************   
   real*8 function Transxa(i,j,k,m)
   integer :: Index , N_Index 
   integer :: i,j,k
   Character*1 :: m   
   real*8 :: Dip_cor, delax_avg, denomenator, press_term 
   real*8 :: Press
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i
      N_Index = Index + 1
      delax_avg = (Deltax(i+1) + Deltax(i)) / 2.0
      Dip_cor = delax_avg ** 2 / ((Zblock(Index)-Zblock(Index+1)) ** 2 + delax_avg ** 2)
      denomenator = 1 / 2.0 *( Deltax(i+1)/(kx(Index+1)* Deltay(j)*Deltaz(Index+1)) + Deltax(i)/(kx(Index)*Deltay(j)*Deltaz(Index)) )
      
      if ( D_Poten(Index , Index + 1 , m ) > 0 ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transxa = betac * Dip_cor * press_term /denomenator      
   
   endfunction Transxa
!***************************************   
   real*8 function Transyb(i,j,k,m) ! grid block (i,j-1,k)
   integer :: Index , N_Index  
   integer :: i,j,k
   Character*1 :: m   
   real*8 :: Dip_cor, delay_avg, denomenator, press_term 
   real*8 :: Press
   
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i
      N_Index = Index - Nx
      delay_avg = (Deltay(j-1) + Deltay(j)) / 2.0
      Dip_cor = delay_avg ** 2 / ( ( Zblock(Index)-Zblock(N_Index) )**2 + delay_avg ** 2)
      denomenator = 1 / 2.0 *( Deltay(j-1)/(ky(N_Index)* Deltax(i)*Deltaz(N_Index)) + Deltay(j)/(ky(Index)*Deltax(i)*Deltaz(Index)) )
      
      if ( D_Poten(Index , N_Index , m ) > 0 ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transyb = betac * Dip_cor * press_term /denomenator   
   endfunction Transyb
!***************************************   
   real*8 function Transya(i,j,k,m) ! grid (i,j+1,k)
   integer :: Index , N_Index 
   integer :: i,j,k
   Character*1 :: m   
   real*8 :: Dip_cor, delay_avg, denomenator, press_term
   real*8 :: Press
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i
      N_Index  = Index + Nx
      delay_avg = (Deltay(j+1) + Deltay(j))/2.0
      Dip_cor = delay_avg ** 2 / ((Zblock(Index)-Zblock( N_Index )) ** 2 + delay_avg ** 2)
      denomenator = 1 / 2.0 *( Deltay(j+1)/(ky( N_Index )* Deltax(i)*Deltaz( N_Index )) + Deltay(j)/(ky(Index)*Deltax(i)*Deltaz(Index)) )
      if ( D_Poten(Index , N_Index , m ) > 0 ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transya = betac * Dip_cor * press_term /denomenator      
   endfunction Transya
!***************************************   
   real*8 function Transzb(i,j,k,m) ! gridblock (i,j,k-1)
   integer :: Index ,  N_Index 
   integer :: i,j,k
   real*8 :: delaz_avg , Press
   Character*1 :: m   
   real*8 :: denomenator, press_term
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i
      N_Index  = Index - Nx * Ny
      denomenator = 1 / 2.0 *( Deltaz( N_Index )/(kz( N_Index )* Deltax(i)*Deltay(j)) + Deltaz(Index)/(kz(Index) * Deltax(i) * Deltay(j)) )
      if ( D_Poten(Index , N_Index , m ) > 0 ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transzb = betac * press_term /denomenator         
   endfunction Transzb         
!***************************************   
   real*8 function Transza(i,j,k,m) ! gridblock (i,j,k+1)
   integer :: Index , N_Index 
   integer :: i,j,k
   real*8 :: delaz_avg , Press
   Character*1 :: m   
   real*8 :: denomenator, press_term
      Index = (k-1) * Nx * Ny + (j-1) * Nx + i
      N_Index  = Index + Nx * Ny
      denomenator = 1 / 2.0 *( Deltaz( N_Index )/(kz( N_Index )* Deltax(i)*Deltay(j)) + Deltaz(Index)/(kz(Index)*Deltax(i)*Deltay(j)) )
      if ( D_Poten(Index , N_Index , m ) > 0 ) then ! upstream block
         if ( m == 'w') then
			press_term = k_rw(Index) / ( mu_w ( Index) * B_w(Index) )            
         Else
			press_term = k_ro(Index) / ( mu_o ( Index) * B_o(Index) )                        
         EndIf     
      else
         if (m == 'w') then
			press_term = k_rw(N_Index) / ( mu_w ( N_Index) * B_w( N_Index ) )                        
         Else
			press_term = k_ro(N_Index) / ( mu_o ( N_Index) * B_o( N_Index ) )                                    
         EndIf           
      Endif 
      Transza = betac * press_term /denomenator         
   endfunction Transza         
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   real*8 function Potential(iii,m)
   
   Character*1 :: m   
   integer:: iii
      if ( m == 'o') then
          Potential = pressure(iii) - density(pressure(iii),m)/144.0 * (Zblock(iii))
      Else
          Potential = pressure(iii) - pco(iii) - density(pressure(iii),m)/144.0 * (Zblock(iii))
      Endif
      
   endfunction Potential            
!##################################################
   subroutine initialize
   ! sets initial pressures according to depth
   Integer :: m ,n , alfa  
   integer :: i,j,k , i6 , iWell, alfawell , ik , B_index
   Real :: Denom , term , Z


If (USEISW) then
	Sw = SWI 
ELSE              
    Do k = 1,Nz
       Do j=1,Ny
          Do i=1,Nx
             alfa = (k-1) * Nx * Ny + (j-1) * Nx + i
             n = size_table + 1
             m = 1
             z = Zblock(alfa)
             do while ( abs(n-m) > 1 )
                i6 = floor(m + (n - m + 1)/2.0) 
                if (z < Z_table(i6) )  then
                    n =  i6 !floor(m2+(n2-m2+1)/2.0);
                else
                    m =  i6 !floor(m2+(n2-m2+1)/2.0);
                end if
             end do
             SW( alfa ) = Sw_table( m ) + ( Zblock(alfa) -  Z_table(m) ) / (  Z_table(n) - Z_table(m) ) * ( Sw_table(n) - Sw_table(m) )                        
          Enddo
       Enddo
    Enddo 	
EndIf


if (USEIP == 'YES') then
		PRESSURE = Initial_P
ELSE
              
    Do k = 1,Nz
       Do j=1,Ny
          Do i=1,Nx
             alfa = (k-1) * Nx * Ny + (j-1) * Nx + i
             n = size_table
             m = 1
             z = Zblock(alfa)
             do while ( abs(n-m) > 1 )
                i6 = floor(m + (n - m + 1)/2.0) 
                if (z < Z_table(i6) )  then
                    n =  i6 !floor(m2+(n2-m2+1)/2.0);
                else
                    m =  i6 !floor(m2+(n2-m2+1)/2.0);
                end if
             end do
             pressure( alfa ) = Op_table( m ) + ( Zblock(alfa) -  Z_table(m) ) / (  Z_table(n) - Z_table(m) ) * ( Op_table(n) - Op_table(m) )  
!!             SW( alfa ) = Sw_table( m ) + ( Zblock(alfa) -  Z_table(m) ) / (  Z_table(n) - Z_table(m) ) * ( Sw_table(n) - Sw_table(m) ) 
                    
          Enddo
       Enddo
    Enddo 
ENDIF    
    
      ! putting pressure and Saturation values into YY vector
    Do k=1,Nz
       Do j=1,Ny
          Do i=1,Nx  
              B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
              alfa = 2 * (B_index ) - 1                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
              YY(alfa) = pressure(B_index)
              YY(alfa + 1) = Sw(B_index)
              
              W_press(B_index) = Pressure(B_index) - PCO(B_index)
              
              k_ro(B_index) = k_r( B_index , 'o')
			  k_rw(B_index) = k_r( B_index , 'w')
			  
              B_o(B_index) = B_f( pressure(B_index) , 'o')
			  B_w(B_index) = B_f( W_press(B_index) , 'w') 	              			  

    		  mu_o(B_index) = mu_f ( pressure(B_index),'o')				
			  mu_w(B_index) = mu_f ( W_press(B_index),'w')				
				
		  
			  
          EndDO
       EndDO
    EndDO
!                Open(1,FILE='initial_pressure.DAT', status='unknown')              
!                !Write(1,*) 'ptable values'
!                Do i = 1 , Ntotal 
!                   Write(1,*) pressure(i) , Zblock(i)
!                EndDO
!                close(1)

            Do i = 1 , Ntotal
				W_press(i) = Pressure(i) - PCO(i)     
            EndDO

   if (.NOT. USEIP == 'YES') then				
       DeAllocate( OP_table )  ! Oil pressure
       DeAllocate( WP_table ) ! Water pressure
       DeAllocate( SW_table )       
       DeAllocate( Z_table) 
   EndIF
       
   Endsubroutine initialize 
!##################################################  

subroutine initialize_BHP
   ! sets initial BHP ; Not too precise but Calc_Delpwf will take care of it
   Integer :: alfa2  , alfa
   Integer :: i , j , i6 , iWell, alfawell  , B_index
   Real*8 :: Denom , term, pres   
   integer :: i_well, i_k,  alfa_well , alfa_w , W_index , N_index
   !Real, Allocatable, Dimension(:) :: O_Volume, W_Volume , mix_den ! O_Volume is the vector of Volume of oil in each completion ( for dt = 1 day , and dp = 1 psi )
   integer :: con 
    
    
    Do i_well = 1 , Nwell       
        W_index = Ntotal + i_well  ! Pwf position in pressure vector
        if ( W_type(i_well) == 'prod') then
           if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! oil rate is given     
              !Do while (con == 0)  
                  term = 0.0
                  Denom = 0.0  
                  Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                 
                     B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                     Denom = Denom + W_index_o(i_well,i_k)
                     term = term + W_index_o(i_well,i_k) * pressure(B_index)
                  EndDo
              !EndDo
              pressure(Ntotal + i_well) = ( term - Schedule(i_well,3,1) )/Denom
		   Elseif ( Schedule(i_well , 2 , sche_N) .EQ. 3 ) then ! Total Liquid rate is given                   
                  term = 0.0
                  Denom = 0.0  
                  Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                 
                     B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                     Denom = Denom + W_index_o(i_well,i_k)+ W_index_w(i_well,i_k)
                     term = term + (W_index_o(i_well,i_k) + W_index_w(i_well,i_k) )* pressure(B_index)
                  EndDo
              pressure(Ntotal + i_well) = ( term - Schedule(i_well,3,1) )/Denom		   
           Else ! BHP was specified
              pressure(Ntotal + i_well) = Schedule(i_well,3,1)
           EndIf
        Else ! injection Well          
           if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! water injection rate is given     
                 term = 0.0
                 Denom = 0.0            
              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                 
                 B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                 Denom = Denom + W_index_w(i_well,i_k)
                 term = term + W_index_w(i_well,i_k) * pressure(B_index)
              EndDo
              pressure(Ntotal + i_well) = ( term + Schedule(i_well,3,1) )/Denom
           Else ! constant BHP injection
              pressure(Ntotal + i_well) = Schedule(i_well,3,1)
           EndIf
        EndIf
        If (pressure(Ntotal + i_well) .LT. MINBHP) then
			pressure(Ntotal + i_well) = MINBHP
		ElseIf (pressure(Ntotal + i_well) .GT. MAXBHP) then
			pressure(Ntotal + i_well) = MAXBHP
		EndIf
    EndDo
   
   Do i=1, Nwell
      B_index = Ntotal + i ! B_index is the block number, while 
      alfa = 2 * Ntotal + i                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
      YY(alfa) = pressure(B_index) 
   EndDo


!                Open(1,FILE='initial_BHP.DAT', status='unknown')              
!                !Write(1,*) 'ptable values'
!                Write(1,*) 'Well Pressure , 1st Block Pressure,    Oil Rate, '
!                Do i_well = 1 , Nwell 
!                   B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( Well_k(i_well , 1) - 1 ) * Nx * Ny
!                   Write(1,*) pressure(i_well) , Schedule(i_well,3,1)
!                EndDO
!                close(1)
   Endsubroutine initialize_BHP
   
!##################################################      
   real*8 function mu_f(p1,m)
   real*8 :: p1
   Integer :: m1,n1 , i5
   Character*1 :: m   
   
   if (m=='w') then
      mu_f = PVTW(4)
   else
        n1 = N_mup
        m1 = 1
        do while (abs(n1-m1).NE.1)
            i5 = floor(m1 + (n1-m1+1)/2.0)
            if (p1 < mu_p(i5,1) )  then
                n1 = floor(m1+(n1-m1+1)/2.0);
            else
                m1 = floor(m1+(n1-m1+1)/2.0);
            end if
        end do
        mu_f = mu_p(m1,2)+ (p1 -mu_p(m1,1))/(mu_p(n1,1)-mu_p(m1,1))*(mu_p(n1,2)-mu_p(m1,2));
   Endif

   Endfunction mu_f
!##################################################     
   real*8 function B_f(p2,m)
   real*8 :: p2 , X
   Integer :: m2,n2 , i4
   Character*1 :: m   
       
       if (m=='o') then
            n2 = N_Bo
            m2 = 1
!           m2 = N_Bo - 1
!           Do
!              if ( p2 > Bo_p(m2,1) )  then
!                 n2 = m2 + 1
!                 exit
!              Else
!                 m2 = m2 - 1
!              EndIf
!           EndDO    
            
            
            do while (abs(n2-m2).NE. 1)
                i4 = FLOOR(m2 + (n2-m2+1)/2.0)
                if (p2 < Bo_p(i4,1) )  then
                    n2 = FLOOR(m2+(n2-m2+1)/2.0);
                else
                    m2 = FLOOR(m2+(n2-m2+1)/2.0);
                end if
            end do
            B_f = Bo_p(m2,2)+ (p2 -Bo_p(m2,1))/(Bo_p(n2,1)-Bo_p(m2,1))*(Bo_p(n2,2)-Bo_p(m2,2));
       Else
            X =  PVTW(3) * (p2 - PVTW(1) )
            B_f = PVTW(2) / ( 1 + X + ( X ** 2 )/2.0 )
       EndIf
       
   Endfunction B_f
   !##################################################   
    real function density(p3,m)  
    
    real*8 :: p3
    Character*1 :: m       
    
        if (m=='o') then
            density = den_o / B_f(p3,'o')
        Else
            density = den_w / B_f(p3,'w')
        Endif
        ! WRITE(*,*) 'DENSITY=', density
        ! pause
    Endfunction density
    
    !##################################################   
  subroutine ref_press_table
 
    implicit none
    
    real :: delta1 
    real*8 :: c1 , c2 , CW1
    real :: maxdepth , mindepth , Delta
    integer :: j1,i1,ref_entry , up_length , down_length ,i , WOC_entry, j
   ! finding maximum depth
       maxdepth = 1.0
       do i1 = 1 , Ntotal 
          !maxdepth = max( maxdepth , Zblock(i1) )
          maxdepth = max( maxdepth , TOPS(i1) + deltaz(i1) )
       End do 
                                                        
        mindepth = 100000
        do i1 = 1 , Ntotal 
            mindepth = min( mindepth , TOPS(i1))
        End do 
                                                        write (*,*) 'minimum depth = ',mindepth
                                                        write (*,*) 'maximum depth = ',maxdepth
                                                        write (*,*) 'ref depth was = ' , initialp(1)

       Delta = 0.05
       
       size_table = FLOOR((maxdepth - mindepth) / Delta) + 100
       
       
       Allocate( OP_table(size_table) )  ! Oil pressure
       Allocate( WP_table(size_table) ) ! Water pressure
       Allocate( SW_table(size_table) )       
       Allocate( Z_table(size_table) )        
       
       ref_entry = abs( FLOOR(( Datum_D - mindepth ) / Delta) ) + 50
       if ( Datum_D < OWC_depth ) then
          Write(*,*) 'The procedure is ok as , Datum_D =',Datum_D,' < ',OWC_depth , '= OWC_depth'
          !pause
       Else
          Write(*,*) 'Sorry, you should change the subroutine for Initial Table'
          write(*,*) 'Datum_D =', Datum_D
          write(*,*)  'OWC_depth=',OWC_depth
          !pause
       EndIf				
       
       Z_table(ref_entry) = Datum_D
       Op_table(ref_entry) = Datum_P ! I supposed that the pressure here is Oil pressure, or in other words the Ref. Press is given above WOC 
       

       up_length = size_table - ref_entry  
       
       Do j1 = 1 , up_length
          Z_table(ref_entry + j1) = Z_table(ref_entry + j1 - 1 ) + Delta 
          if  ( Z_table( ref_entry + j1 ) > OWC_depth ) then 
             Z_table( ref_entry + j1 ) =  OWC_depth
             Delta = OWC_depth - z_table(ref_entry + j1 - 1 )
             c1 = Op_table( ref_entry + j1 - 1) + density( Op_table(ref_entry + j1 - 1) , 'o' ) / 144.0 *  Delta 
             delta1 = 1d-1
             do while ( delta1 > 1d-4 )
                c2 = Op_table( ref_entry + j1 - 1) + ( density( Op_table(ref_entry + j1 - 1) , 'o' ) + density( c1, 'o' ) ) /2.0  / 144.0 * Delta 
                delta1 = abs(c2 - c1)
                c1 = c2
             end do
             WOC_entry = ref_entry  + j1
             Op_table( ref_entry  + j1) = c1
             Wp_table( ref_entry  + j1) = Op_table( ref_entry  + j1) - OWC_PCOW


             Exit
          EndIF
          
          c1 = Op_table( ref_entry + j1 - 1) + density( Op_table(ref_entry + j1 - 1) , 'o' ) / 144.0 *  Delta 
          !write(*,*) density( p_table(ref_entry + j1-1,2) )
          !write(*,*) p_table(ref_entry + j1-1,2)          
                    
          delta1 = 1d-1
          do while ( delta1 > 1d-6 )
             c2 = Op_table(ref_entry + j1 - 1) + ( density( Op_table(ref_entry + j1 - 1) , 'o' ) + density( c1, 'o' ) ) /2.0 / 144.0 * Delta !( p_table(ref_entry + j1,1) - p_table(ref_entry + j1-1,1) )     density( ( Op_table(ref_entry + j1 - 1) + c1 ) / 2.0 , 'o') 
             delta1 = abs(c2 - c1)
             !write(*,*) 'p_table =' ,              
             c1 = c2
          end do
          Op_table( ref_entry  + j1) = c1
          
!                                    write (*,*) Op_table( ref_entry  + j1)
!                                    pause          
       End DO
       
       
       Delta = 0.05
	If (.NOT. USEISW) then       
		   Do i = WOC_entry + 1 , size_table !up_length  ! below OWC there is only water , which means Oil pressure should be calculated according to water pressure calculated first
			  Z_table(i) = Z_table( i - 1 ) + Delta 
			  CW1 = Wp_table( i - 1) + density( Wp_table(i - 1) , 'w' ) / 144.0 *  Delta
	                              
			  delta1 = 1d-1
			  do while ( delta1 > 1d-6 )
				 c2 = Wp_table(i - 1) + ( density( Wp_table(i - 1) , 'w') + density( CW1 , 'w' ) )/2.0 / 144.0 * Delta !( p_table(ref_entry + j1,1) - p_table(ref_entry + j1-1,1) )
				 delta1 = abs(c2 - CW1)
				 !write(*,*) 'p_table =' ,              
				 CW1 = c2
			  end do
			  Wp_table(i) = CW1
			  Op_table(i) = Wp_table(i) +  OWC_PCOW
		   EndDO
	EndIF
       Delta = 0.05  
       down_length = ref_entry - 1
       
       Do j1 = 1 , down_length
          i = ref_entry - j1
          Z_table(i) = Z_table( i + 1 ) - Delta 
          c1 = Op_table(i + 1) - density( Op_table(i + 1 ) , 'o' ) / 144.0 * Delta !( p_table(ref_entry - j1,1) - p_table(ref_entry - j1 + 1,1) )
          delta1 = 1d-1
          do while ( delta1 > 1d-6 )
             c2 = Op_table( i + 1) - ( density( Op_table(i + 1) , 'o' ) + density( c1 , 'o' ) ) / 2.0 / 144.0 * Delta !( p_table(ref_entry - j1,1) - p_table(ref_entry - j1 + 1,1) )
             delta1 = abs(c2 - c1)
             c1 = c2
          enddo
          Op_table(i) = c1          
       End DO
       
       ! PW_table values before WOC are calculated at this step
		If (WOC_entry >	size_table) then
				WOC_entry =	size_table
		EndIf
       Do i = 1 , WOC_entry - 1 
          j = WOC_entry - i
          Delta = Z_table(j + 1) - Z_table( j ) 
          CW1 = Wp_table( j + 1) - density( Wp_table(j + 1) , 'w' ) / 144.0 *  Delta
                              
          delta1 = 1d-1
          do while ( delta1 > 1d-4 )
             c2 = Wp_table(j + 1) - ( density( Wp_table(j + 1) , 'w' ) + density( CW1 , 'w' ) ) / 2.0 / 144.0 * Delta !( p_table(ref_entry + j1,1) - p_table(ref_entry + j1-1,1) )
             delta1 = abs(c2 - CW1)
             CW1 = c2
          end do
          Wp_table(j) = CW1       
       EndDo
              
       !write(*,*) down_length 
       !write(*,*) size_p_table - up_length - 1
       !pause
       
              
                Open(1,FILE='OP_table.DAT', status='unknown')              
                Write(1,*) 'Depth ,    Oil Pressure,  Water Pressure,          Sw'
                
                Do j1 = 1, size_table
                   Write(1,*)  Z_table(j1) , Op_table(j1) , Wp_table(j1)
                EndDo   !                
                close(1)
                
       If (.NOT. USEISW) then
		   Do i = 1 , size_table!WOC_entry
			   Sw_table ( i ) = PCO_inv( Op_table(i) - Wp_table(i) ) 
		   EndDO
	   EndIf
!       write(*,*) 'maximum p_table is ' , max(op_table),'   min P_table is ', min(op_table)
!       pause
              
!                Open(1,FILE='Sw_table.DAT', status='unknown')              
!                Write(1,*) 'Depth ,    Oil Pressure - Water Pressure,          Sw'
!                
!                Do j1 = 1, size_table
!                   Write(1,*)  Z_table(j1) , Op_table(j1) - Wp_table(j1), Sw_table(j1)
!                EndDo   
!                
!                close(1)

       
    Endsubroutine ref_press_table
   !##################################################   
    real function c_f(p4)  
    real*8 :: p4
        c_f = 0.00005 / B_f(p4,'o')
    Endfunction c_f
    
    !##################################################       
    real*8 function PCO(ind) ! computed the cappilary pressure in the given block according to its saturation
    
    real :: swater
    Integer :: m,n , ind
            swater = sw(ind)
            n = n_Kr
            m=1
            do while (abs(n-m) .NE. 1)
                if (swater < Sw_mat( floor( m+(n-m+1)/2.0 ) ))  then
                    n = floor(m+ ( n - m + 1 ) / 2.0 );
                else
                    m = floor(m + ( n - m + 1) / 2.0);
                end if
            end do
            PCO = Pcow(m)+ (swater - Sw_mat(m))/(Sw_mat(n) - Sw_mat(m))*(Pcow(n) - Pcow(m));
    !PCO = - 34.98 * Swater ** 3 + 83.96 * Swater**2 - 67.17 * Swater
    Endfunction PCO


    !##################################################       
        
    real function k_r(ind , m1 )
    ! You Must Notice that in this Function the assumption is that the value of kro is 0 for the last 2 values of Sw is Sw_mat matrix
    ! E.G. Sw = 0.8 , kro = 0 ; Sw = 1 , kro =0 , 
    ! if kro(Sw = Sw_mat(n_kr - 1 ) ) .NE. 0 , then this function should be revised    

    !-------------------------------------Variable Declaration Part-------------------------------------------------------
    Implicit None
    !
    Character*1 :: m1  ! 'oil or water'
    real :: swater
    Integer :: m , n , ind
            swater = sw(ind)
            n = n_Kr
            m = 1

!		   m = n_Kr - 1
!           Do
!              if ( swater > Sw_mat(m) )  then   ! It means that Water Saturaion is never below ( even by a very very small value) the minimum value
!                 n = m + 1
!                 exit
!              Else
!                 m = m - 1
!              EndIf
!           EndDO    

            do while (abs(n-m) .NE. 1)
                if (swater < Sw_mat(floor( m + (n-m+1)/2.0) ) )  then
                    n = floor(m + ( n - m + 1) / 2.0);
                else
                    m = floor(m + ( n - m + 1) / 2.0);
                end if
            end do
            if  (m1 == 'o') then
				if (swater .LE. Sw_mat(1) ) then
					k_r = kro(1)
				Elseif (swater .GE. Sw_mat(n_Kr) ) then             
					k_r = 0.0
				Else
					k_r = kro(m)+ (swater -Sw_mat(m))/(Sw_mat(n)-Sw_mat(m))*(kro(n)-kro(m))
				EndIf
            else
				if (swater .LE. Sw_mat(1)) then
!					sw(ind) = Sw_mat(1)
					k_r = 0.0
				Elseif (swater .GE. 1.0 ) then             			
					k_r = krw(n_Kr)
				Else            
					k_r = krw(m)+ (swater -Sw_mat(m))/(Sw_mat(n)-Sw_mat(m))*(krw(n)-krw(m));
				EndIf
            EndIf
    endfunction k_r
    !##################################################       
        
!    real*8 function W_press(ind) ! this function calculates Water pressure in Block # ind by subtracting Capillary pressure from Oil Pressure
!    Implicit none
!    integer :: ind
!    real :: Swater
!    
!        Swater = SW(ind)
!        W_press = Pressure(ind) - PCO(ind)
!    Endfunction W_press 
    !##################################################       
  real function PCO_inv(S)
    
    real*8 :: S ! the type of the input variable is PRESSURE 
    Integer :: m,n , ind
        n = n_Kr
        m=1
        if ((S .LE. Pcow(n_kr))  ) then ! minimum Capillary pressure  .OR. (abs( S - Pcow(n_kr)) .LE. 0.00004 )
           PCO_inv = Sw_mat(n_kr) !- Epsilon
        Elseif ( S .GE. Pcow(1) ) then ! maximum capillary pressure
           PCO_inv = Sw_mat(1) 
        Else
           m = n_kr - 1
           Do
              if ( S < Pcow(m) )  then
                 n = m + 1
                 exit
              Else
                 m = m - 1
              EndIf
           EndDO    
           PCO_inv = Sw_mat(m) + ( S - Pcow(m) ) / ( Pcow(n) - Pcow(m) ) * ( Sw_mat(n) - Sw_mat(m) )
        EndIf
                    
    Endfunction PCO_inv
    !##################################################       

   subroutine initial_swi
   ! sets initial pressures according to depth
   Integer :: m ,n , alfa  
   integer :: i,j,k , i6 , iWell, alfawell , ik , B_index
   Real :: Denom , term , Z
   !Allocate(pressure(NTOTAL + Nwell ))
       
    Do k = 1,Nz
       Do j=1,Ny
          Do i=1,Nx
             alfa = (k-1) * Nx * Ny + (j-1) * Nx + i
             n = size_table + 1
             m = 1
             z = Zblock(alfa)
             do while ( abs(n-m) > 1 )
                i6 = floor(m + (n - m + 1)/2.0) 
                if (z < Z_table(i6) )  then
                    n =  i6 !floor(m2+(n2-m2+1)/2.0);
                else
                    m =  i6 !floor(m2+(n2-m2+1)/2.0);
                end if
             end do
             pressure( alfa ) = Op_table( m ) + ( Zblock(alfa) -  Z_table(m) ) / (  Z_table(n) - Z_table(m) ) * ( Op_table(n) - Op_table(m) )  
             SW( alfa ) = 0.2 
     
          Enddo
       Enddo
    Enddo 
      ! putting pressure and Saturation values into YY vector
    Do k=1,Nz
       Do j=1,Ny
          Do i=1,Nx  
              B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
              alfa = 2 * (B_index ) - 1                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
              YY(alfa) = pressure(B_index)
              YY(alfa + 1) = Sw(B_index)
          EndDO
       EndDO
    EndDO
                Open(1,FILE='initial_pressure.DAT', status='unknown')              
                !Write(1,*) 'ptable values'
                Do i = 1 , Ntotal 
                   Write(1,*) pressure(i) , Zblock(i)
                EndDO
                close(1)

   Endsubroutine initial_swi 
!##################################################  
   real*8 function D_Poten(ind1,ind2,m)
   
   Character*1 :: m   
   integer:: ind1, ind2
   real*8 :: den_avg
   
      if ( m == 'o') then
		 den_avg = ( density( pressure(ind1),'o' ) + density( pressure(ind2),'o' ) )/ 2.0
         D_Poten = (pressure(ind1) - pressure(ind2)) - den_avg / 144.0 * ( Zblock(ind1) - Zblock(ind2) )
      Else
		 den_avg = ( density( W_press(ind1) ,'w') + density( W_press(ind2) , 'w') )/ 2.0
         D_Poten = W_press(ind1) - W_press(ind2) - den_avg / 144.0 * ( Zblock(ind1) - Zblock(ind2) )
      Endif
      
   endfunction D_Poten           
!*******************************************************************************************
   subroutine change_Sw
   
   integer :: i , j , k , B_index, counter  
	  
	  counter = 0 
      Do k = 1 , Nz
         Do j = 1 , Ny
            Do i = 1 , Nx          
                B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
                if ( abs(1 - Sw(B_index)) < 1d-4) then
                   Sw(B_index) = 1 - 1d-4
                   counter = counter + 1
                EndIf
            EndDO
         EndDO
      EndDO   
   
   Endsubroutine change_Sw
!*******************************************************************************************   
    !##################################################       
    real function PCO_Sw(W_saturation) ! computed the cappilary pressure in the given block according to its saturation
    
    real*8 :: W_saturation
    Integer :: m,n 

            n = n_Kr
            m=1
            do while (abs(n-m)>1)
                if (W_saturation < Sw_mat( floor( m+(n-m+1)/2.0 ) ))  then
                    n = floor(m+ ( n-m+1 ) / 2.0 );
                else
                    m = floor(m + ( n - m + 1) / 2.0);
                end if
            end do
            PCO_Sw = Pcow(m)+ (W_saturation - Sw_mat(m))/(Sw_mat(n) - Sw_mat(m))*(Pcow(n) - Pcow(m));
        
    Endfunction PCO_Sw


    !##################################################               
     
EndModule Trans