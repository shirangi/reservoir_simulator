
module Mat_Calc

!!
Use globalvar
Use Trans
Use well
Use  Derivative

! in this Module RHS and Jacobean matrix are defined
contains 

Subroutine MATRIX



	call calc_RHS
	call calc_Jac


EndSubroutine MATRIX




Subroutine calc_RHS    
!-------------------------------------Variable Declaration Part-------------------------------------------------------
Implicit None

integer:: alfa , alfa_N, i ,j ,k , i_well , B_index , N_index , W_index , alfa_WB , alfa_W	, ii, jj, kk

real*8 :: Vb , den_avg ,  CC , chek , chek2, chek3 , X
real*8 ::Txbo , Txao , Tybo , Tyao , Tzbo , Tzao  ,Txbw , Txaw , Tybw , Tyaw , Tzbw , Tzaw
real*8 :: acc_w, acc_o , poros_n , poros_p , acc_term , Well_term , pr , pr2 , Dxpb , Dxsb
integer :: i_k , alfa_wprim , alfa_wpp , ycounter ,xcounter , WB_index, B_alfa
!-------------------------------------*************************-------------------------------------------------------

	FCN = 0	
                                                                                                  
	! Inner LOOP
	Do k= 1 , Nz 
	   Do j= 1 , Ny 
		  Do i= 1 , Nx	
			 B_index = (k-1) * Nx * Ny + (j-1) * Nx + i
			 alfa = 2 * ( B_index -1 ) +1

			 Txbo = 0
			 Txao = 0
			 Tybo = 0
			 Tyao = 0
			 Tzbo = 0
			 Tzao = 0
			 Txbw = 0
			 Txaw = 0
			 Tybw = 0
			 Tyaw = 0
			 Tzbw = 0
			 Tzaw = 0
	         

!			^^^^^^^^^^^^^^^^^^^^^^^^^^^
			 if (i .NE. 1) then
				Txbo = Transxb(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Txbo *  D_Poten(B_index - 1 , B_index , 'o') 
				Txbw = Transxb(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Txbw *  D_Poten(B_index - 1 , B_index , 'w') 
			 end if 			
	         
!			^^^^^^^^^^^^^^^^^^^^^^^^^^^
			 if (i .NE. Nx) then
				Txao = Transxa(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Txao *  D_Poten(B_index + 1 , B_index , 'o') 
				Txaw = Transxa(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Txaw *  D_Poten(B_index + 1 , B_index , 'w') 
			 end if 			

!			^^^^^^^^^^^^^^^^^^^^^^^^^^^	         
			 if (j .NE. 1) then
				Tybo = Transyb(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Tybo *  D_Poten(B_index - Nx , B_index , 'o') 
				Tybw = Transyb(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Tybw * D_Poten(B_index - Nx , B_index , 'w') 
			 end if 
	         
!			^^^^^^^^^^^^^^^^^^^^^^^^^^^
			 if (j .NE. Ny) then
				Tyao = Transya(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Tyao * D_Poten(B_index + Nx , B_index , 'o') 
				Tyaw = Transya(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Tyaw * D_Poten(B_index + Nx , B_index , 'w') 
			 end if 
			
!			^^^^^^^^^^^^^^^^^^^^^^^^^^^	         
			 if (k .NE. 1) then
				N_index = B_index - Nx * Ny
				Tzbo = Transzb(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Tzbo * D_Poten( N_index , B_index , 'o') 
				Tzbw = Transzb(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Tzbw * D_Poten( N_index , B_index , 'w')  
			 end if 

!			^^^^^^^^^^^^^^^^^^^^^^^^^^^
			 if (k .NE. Nz) then
				N_index = B_index + Nx * Ny
				Tzao = Transza(i,j,k,'o') 
				FCN(alfa) = FCN(alfa) + Tzao * D_Poten( N_index , B_index , 'o')  
				Tzaw = Transza(i,j,k,'w') 
				FCN( alfa + 1 ) = FCN( alfa + 1 ) + Tzaw * D_Poten( N_index , B_index , 'w') 
			 end if 

!			^^^^^^^^^^^^^^^^^^^^^^^^^^^	         
			 !C_ROCK , P_ref 
			 X = C_Rock * ( pressure(B_index ) - P_ref )
			 poros_n  = phi(B_index) * (1 + X + 1.0 / 2.0 * X**2  ) 
			 
			 X = C_Rock * ( past_p(B_index) - P_ref )
			 poros_p  = phi(B_index) * (1 + X + 1.0 / 2.0 * X**2  ) 
			 
			 Vb = Deltax(i) * Deltay(j)* Deltaz(B_index)
	         
	         
			 acc_o = Vb * poros_n  * (1 - Sw(B_index)) / ( alfac * B_o(B_index) * d_time) - Vb * poros_p * (1 - past_Sw(B_index)) / ( alfac * B_f( past_p(B_index),'o') * d_time)
			 acc_w = Vb * poros_n  *  Sw(B_index) / ( alfac * B_w(B_index) * d_time) - Vb * poros_p  * past_Sw(B_index) / ( alfac * B_f(past_pw(B_index),'w') * d_time)
			 
			 FCN(alfa) = FCN(alfa) - acc_o
			 FCN(alfa + 1) = FCN(alfa + 1) - acc_w

		  Enddo
	   Enddo
	Enddo
	

	 ! I supposed that the constrants can just be Oil rate or BHP in case of a production well, so There is no case for constant Total Liquid production
	Do i_well = 1 , Nwell
		W_index = Ntotal + i_well  ! Pwf position in pressure vector
		alfa_w = 2 * Ntotal + i_well 
		if ( W_type(i_well) == 'prod') then
		   if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! oil rate is given
			   FCN(alfa_w) = FCN(alfa_w) - Schedule(i_well , 3 ,sche_N) ! Well equation 
			   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
				   WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2 ) - 1 ) * Nx + (i_k - 1) * Nx * Ny ! Well block index
				   alfa_WB = 2 * WB_index - 1   ! YY index of the well block
				   ! f(o)
				   FCN( alfa_WB ) = FCN( alfa_WB ) - W_index_o(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				   ! f(w)
				   FCN( alfa_WB + 1 ) = FCN( alfa_WB + 1 ) - W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )                   

				   ! Fwell
				   FCN(alfa_w) = FCN(alfa_w) + W_index_o(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )     				   
			   EndDO       
		   ELSEif ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then ! BHP is given
			   FCN(alfa_w) = 0.0 ! Well equation 
			   YY (alfa_w) = Schedule(i_well , 3 ,sche_N)
			   pressure(W_index) = Schedule(i_well , 3 ,sche_N)
			   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
				   WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2 ) - 1 ) * Nx + (i_k - 1) * Nx * Ny ! Well block index
				   alfa_WB = 2 * WB_index - 1   ! YY index of the well block
				   ! f(o)
				   FCN( alfa_WB ) = FCN( alfa_WB ) - W_index_o(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				   ! f(w)
				   FCN( alfa_WB + 1 ) = FCN( alfa_WB + 1 ) - W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )    
			   EndDO       
		   Elseif ( Schedule(i_well , 2 , sche_N) .EQ. 3 ) then ! Total Liquid Rate Specified
			   FCN(alfa_w) = FCN(alfa_w) - Schedule(i_well , 3 ,sche_N) ! Well equation 
			   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
				   WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2 ) - 1 ) * Nx + (i_k - 1) * Nx * Ny ! Well block index
				   alfa_WB = 2 * WB_index - 1   ! YY index of the well block
				   ! f(o)
				   FCN( alfa_WB ) = FCN( alfa_WB ) - W_index_o(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				   ! f(w)
				   FCN( alfa_WB + 1 ) = FCN( alfa_WB + 1 ) - W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )                   

				   ! Fwell
				   FCN(alfa_w) = FCN(alfa_w) + (W_index_o(i_well,i_k) + W_index_w(i_well,i_k) ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )     				   
			   EndDO       			   
		   EndIf
		Else   ! W_type(i_well) == 'inje'
		   if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! water injection rate is given
			   FCN(alfa_w) = FCN(alfa_w) + Schedule(i_well , 3 ,sche_N) ! Well equation  !!! NOTICE that for injection, the equation is the same but, RATE is NEGATIVE
			   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
				   WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2 ) - 1 ) * Nx + (i_k - 1) * Nx * Ny ! Well block index
				   alfa_WB = 2 * WB_index - 1   ! YY index of the well block
				   ! No Oil is entering the Block
				   ! fwell
				   FCN( alfa_w ) = FCN( alfa_w ) + W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				   ! f(w)
				   FCN( alfa_WB + 1 ) = FCN( alfa_WB + 1 ) - W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )    
			   EndDO       
		   ELSEif ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then ! BHP is gibven
			   FCN(alfa_w) = 0 ! Well equation 
			   YY (alfa_w) = Schedule(i_well , 3 ,sche_N)
			   pressure(W_index) = Schedule(i_well , 3 ,sche_N)
			   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
				   WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2 ) - 1 ) * Nx + (i_k - 1) * Nx * Ny ! Well block index
				   alfa_WB = 2 * WB_index - 1   ! YY index of the well block
				   ! f(w)
				   FCN( alfa_WB + 1 ) = FCN( alfa_WB + 1 ) - W_index_w(i_well,i_k) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )    
			   EndDO       
		   EndIf	               
		Endif
	EndDo

	maxFCN = maxval(abs(FCN))
	maxentry = maxloc(abs(FCN))  
	
!			Open(3,FILE='RHS.DAT', status='unknown')
!			write(3,*) 'index ,  i   ,  j   , Entry Value '
!			Do i=1,N_equ
!				write(3,*) i , FCN(i)
!			EndDo
!			close(3)  

EndSubroutine calc_RHS


! ##########################################################################################################################################################



Subroutine calc_Jac

!-------------------------------------Variable Declaration Part-------------------------------------------------------
Implicit None

integer:: alfa , alfa_N, i ,j ,k , i_well , B_index , N_index , W_index , alfa_WB , alfa_W

real*8 :: Vb , den_avg ,  CC , chek , chek2, chek3 , X
real*8 ::Txbo , Txao , Tybo , Tyao , Tzbo , Tzao  ,Txbw , Txaw , Tybw , Tyaw , Tzbw , Tzaw , Del_PO, Del_PW
real*8 :: acc_w, acc_o , poros_n , poros_p , acc_term , Well_term , pr , pr2 , Dxpb , Dxsb
integer :: i_k , alfa_wprim , alfa_wpp , ycounter ,xcounter , WB_index, B_alfa , Jac_index
Real*8 , Allocatable, Dimension(:) :: Jac , i_Jac , j_Jac !! Instead of Big matrix Jacobian!
real*8 :: term
! Jacobian has at most N * 15 Entries ( N is the number of equations ), because of 6 neighbour blocks , 1 Center term , and 1 well term
! However Many of these entries are zero
! 1st Entry = Po k-1 , 
! 2nd Entry = Sw k-1 , 
! 3rd Entry = Po j-1, 
! 4th Entry = Sw j-1
! 5th Entry = Po i-1 
! 6th Entry = Sw i-1 
! 7th Entry = Po Center 
! 8th Entry = Sw Center
! 9th Entry = Po i + 1
! 10th Entry= Sw i + 1
! 11th Entry= Po j + 1
! 12th Entry= Sw j + 1
! 13th Entry= Po k + 1
! 14th Entry= Sw k + 1
! 15th Entry= Well ( Pwf)
!-------------------------------------*************************-------------------------------------------------------
	ALLOCATE (Jac(N_equ*15) ,i_Jac(N_equ*15) ,j_Jac(N_equ*15)  )	 
	call change_Sw
	Jac = 0 ; i_Jac = 0 ; j_Jac = 0

	! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^  JACOBIAN GENERATION  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


	Do k=1,Nz
	   Do j=1,Ny
		  Do i=1,Nx
	      
			 B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
			 alfa = 2 * (B_index ) - 1                 ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
													  ! in pressure vector
													  ! Since we have the pressure matrix each time ( pressure values extracted from Y vector and put into Pressure vector )			 
			
			 if (k .NE. 1) then 
	             
				 N_index = B_index - Nx * Ny    ! Neighbour block Index
				 alfa_N = alfa - 2 * Nx * Ny     ! Neighbour index in Y vector           
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')				 
				 
				 Tzbo =	Transzb(i,j,k,'o')
				 Tzbw =	Transzb(i,j,k,'w')				 				 				 

				 ! d(fo(ijk))/d(po(ijk-1))  , P(ijk-1) = Y (alfa_N)
				 term = DTrzb_DP(i,j,k,'o','b') * Del_PO + Tzbo * ( 1.0 - gama_der( N_index , 'o' )*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 1
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N		
				 				 
				 ! d(fw(ijk))/d(po(ijk-1))
				 term = DTrzb_DP(i,j,k,'w','b') * Del_PW + Tzbw * ( 1.0 - gama_der( N_index , 'w') * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 1
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N					 

				 ! d(fo(ijk))/d(Sw(ijk-1))  , Sw(ijk-1) = Y (alfa_N + 1)
				 term = DTrzb_DS(i,j,k,'o','b') * Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 2
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa_N + 1					 
				 				 
				 ! d(fw(ijk))/d(Sw(ijk-1))
				 term = DTrzb_DS(i,j,k,'w','b') * Del_PW + Tzbw * PCO_der(N_index) * ( - 1.0 + gama_der(N_index , 'w' ) *( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 2
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1					 
	             
				 ! d(fo(ijk))/d(po(ijk))
				 term =  DTrzb_DP(i,j,k,'o','c') * Del_PO + Tzbo * ( - 1.0 - gama_der(B_index , 'o') * ( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 				 				 
				 ! d(fw(ijk))/d(po(ijk))
				 term =  DTrzb_DP(i,j,k,'w','c') * Del_PW + Tzbw * ( - 1.0 - gama_der( B_index , 'w' ) * ( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa	
				 
				 ! d(fo(ijk))/d(Sw(ijk))
				 term =  DTrzb_DS(i,j,k,'o','c') * Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term =  DTrzb_DS(i,j,k,'w','c') * Del_PW + Tzbw * PCO_der(B_index) *( 1.0 + gama_der( B_index , 'w' ) * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1					 

			 EndIf
			 if (j .NE. 1) then   
				 N_index = B_index - Nx    ! Neighbour block Index
				 alfa_N = alfa - 2 * Nx     ! Neighbour index in Y vector
				 
				 Tybo =	Transyb(i,j,k,'o')
				 Tybw =	Transyb(i,j,k,'w')				 
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')				 

				 ! d(fo(ijk))/d(po(ij-1k))  , P(ij-1k) = Y (alfa_N)
				 term = DTryb_DP(i,j,k,'o','b') *  Del_PO + Tybo * ( 1.0 - gama_der( N_index , 'o' )*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 3
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N	
				 				 
				 ! d(fw(ijk))/d(po(ij-1k))
				 term = DTryb_DP(i,j,k,'w','b') * Del_PW + Tybw * ( 1.0 - gama_der( N_index , 'w')*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 3
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N					 

				 ! d(fo(ijk))/d(Sw(ij-1k))  , Sw(ij-1k) = Y (alfa_N + 1 )
				 term = DTryb_DS(i,j,k,'o','b') * Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 4
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa_N + 1				 
				 
				 ! d(fw(ijk))/d(Sw(ij-1k))
				 term = DTryb_DS(i,j,k,'w','b') * Del_PW +  Tybw * ( - PCO_der(N_index) )*( 1.0 - gama_der(N_index , 'w' ) *( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 4
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1				 
	             
				 ! d(fo(ijk))/d(po(ijk))
				 term = DTryb_DP(i,j,k,'o','c') * Del_PO + Tybo * ( -1.0 - gama_der( B_index , 'o') * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 				 
				 ! d(fw(ijk))/d(po(ijk))
				 term = DTryb_DP(i,j,k,'w','c') * Del_PW + Tybw * ( -1.0 - gama_der( B_index , 'w')*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa					 
	             
				 ! d(fo(ijk))/d(Sw(ijk))
				 term = DTryb_DS(i,j,k,'o','c') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term = DTryb_DS(i,j,k,'w','c') * Del_PW + Tybw * PCO_der(B_index) * ( 1.0 + gama_der( B_index , 'w' ) * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1					 

			 EndIf
			 
			 if (i .NE. 1) then  
				 N_index = B_index - 1 
				 alfa_N = alfa - 2      ! Neighbour index in X vector
				 
				 Txbo =	Transxb(i,j,k,'o')
				 Txbw =	Transxb(i,j,k,'w')				 
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')	
				 			 				 
				 ! d(fo(ijk))/d(po(i-1jk))   , p(i-1) = Y ( i-2)
				 term = DTrxb_DP(i,j,k,'o','b') * Del_PO + Txbo * ( 1.0 - gama_der(B_index - 1 , 'o' ) * ( Zblock(B_index-1) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 5
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N				 
				 
				 ! d(fw(ijk))/d(po(i-1jk))
				 term = DTrxb_DP(i,j,k,'w','b') * Del_PW + Txbw * ( 1.0 - gama_der(B_index - 1 , 'w' )*( Zblock(B_index-1) - Zblock(B_index) ) )
				 Jac_index = 15 * ( alfa + 1 - 1) + 5
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N				 

				 ! d(fo(ijk))/d(Sw(i-1jk))   , Sw(i-1) = Y ( alfa-1)
				 term = DTrxb_DS(i,j,k,'o','b') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 6
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N + 1
				 
				 ! d(fw(ijk))/d(Sw(i-1jk))
				 term = DTrxb_DS(i,j,k,'w','b') * Del_PW + Txbw * PCO_der(N_index) * ( - 1.0 + gama_der(N_index , 'w') * ( Zblock( B_index - 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 6
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1

				 ! d(fo(ijk))/d(po(ijk))
				 term = DTrxb_DP(i,j,k,'o','c') * Del_PO + Txbo * ( -1.0 - gama_der(B_index , 'o' ) * ( Zblock( B_index - 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 
				 ! d(fw(ijk))/d(po(ijk))
				 term = DTrxb_DP(i,j,k,'w','c') * Del_PW + Txbw  *  (-1.0 - gama_der(B_index , 'w' ) * ( Zblock( B_index - 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa				 

				 ! d(fo(ijk))/d(Sw(ijk))
				 term = DTrxb_DS(i,j,k,'o','c') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term =  DTrxb_DS(i,j,k,'w','c') * Del_PW + Txbw  * PCO_der(B_index) * ( 1.0 + gama_der(B_index , 'w' ) * ( Zblock( B_index - 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1					 
				 
	             
			 EndIf

			 if (i .NE. Nx) then   
				 N_index = B_index + 1    
				 alfa_N = alfa + 2
				 Txao =	Transxa(i,j,k,'o')
				 Txaw =	Transxa(i,j,k,'w')				 				 
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')				 				 				 
				 
				 ! d(fo(ijk))/d(po(ijk))
				 term =  DTrxa_DP(i,j,k,'o','c') * Del_PO + Txao * (-1.0 - gama_der( B_index , 'o' ) * ( Zblock( B_index + 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 				 
				 ! d(fw(ijk))/d(po(ijk))
				 term = DTrxa_DP(i,j,k,'w','c') * Del_PW + Txaw * (-1.0 - gama_der( B_index , 'w' ) * ( Zblock( B_index + 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa						 

				 ! d(fo(ijk))/d(Sw(ijk))
				 term = DTrxa_DS(i,j,k,'o','c') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term = DTrxa_DS(i,j,k,'w','c') * Del_PW + Txaw * PCO_der(B_index) * ( 1.0 + gama_der( B_index , 'w' ) * ( Zblock( B_index + 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1					 

				 ! d(fo(ijk))/d(po(i+1jk))   , P(i+1) = Y ( i + 2)
				 term = DTrxa_DP(i,j,k,'o','a') *  Del_PO + Txao * ( 1.0 - gama_der(B_index + 1, 'o' ) * ( Zblock( B_index + 1 ) - Zblock( B_index ) ) ) ! 
				 Jac_index = 15 * ( alfa - 1) + 9
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N				 
				 				 
				 ! d(fw(ijk))/d(po(i+1jk))
				 term = DTrxa_DP(i,j,k,'w','a') *  Del_PW + Txaw * ( 1.0 - gama_der( B_index + 1 , 'w' ) * ( Zblock( B_index + 1 ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 9
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N				 
				 

				 ! d(fo(ijk))/d(Sw(i+1jk))   , Sw(i+1) = Y ( i + 3)
				 term = DTrxa_DS(i,j,k,'o','a') *  Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 10
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa_N + 1
				 
				 ! d(fw(ijk))/d(Sw(i+1jk))
				 term = DTrxa_DS(i,j,k,'w','a') *  Del_PW +  Txaw * PCO_der(N_index) * ( -1.0  + gama_der(N_index , 'w' )*( Zblock( N_index ) - Zblock( B_index ) ) ) ! 
				 Jac_index = 15 * ( alfa + 1 - 1) + 10
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1
	                
			 EndIf


			 if (j .NE. Ny) then   
				 N_index = B_index + Nx    ! Neighbour block Index
				 alfa_N = alfa + 2 * Nx     ! Neighbour index in Y vector
				 
				 Tyao =	Transya(i,j,k,'o')
				 Tyaw =	Transya(i,j,k,'w')				 				 
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')				 

				 ! d(fo(ijk))/d(po(ijk))
				 term = DTrya_DP(i,j,k,'o','c') * Del_PO + Tyao * (-1.0 - gama_der(B_index , 'o' )*( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 				 
				 ! d(fw(ijk))/d(po(ijk))
				 term = DTrya_DP(i,j,k,'w','c') * Del_PW + Tyaw * (-1.0 - gama_der(B_index , 'w')*( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa						 

				 ! d(fo(ijk))/d(Sw(ijk))
				 term = DTrya_DS(i,j,k,'o','c') * Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term = DTrya_DS(i,j,k,'w','c') * Del_PW + Tyaw * PCO_der(B_index) * ( 1.0 + gama_der( B_index , 'w' ) * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1						 

				 ! d(fo(ijk))/d(po(ij+1k))  , P(ij+1k) = Y (alfa_N)
				 term =  DTrya_DP(i,j,k,'o','a') * Del_PO + Tyao * ( 1.0 - gama_der( N_index , 'o' )*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 11
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N				 
				 
				 ! d(fw(ijk))/d(po(ij+1k))
				 term =  DTrya_DP(i,j,k,'w','a') * Del_PW + Tyaw * ( 1.0 - gama_der( N_index , 'w')*( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 11
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N					 

				 ! d(fo(ijk))/d(Sw(ij+1k))  , Sw(ij+1k) = Y (alfa_N + 1)
				 term = DTrya_DS(i,j,k,'o','a') * Del_PO 
				 Jac_index = 15 * ( alfa - 1) + 12
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa_N + 1					 
				 
				 ! d(fw(ijk))/d(Sw(ij+1k))
				 term = DTrya_DS(i,j,k,'w','a') * Del_PW +  Tyaw * PCO_der(N_index) * ( -1.0 + gama_der(N_index , 'w' ) * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 12
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1				 
	             
			 EndIf


			 if (k .NE. Nz) then   
				 N_index = B_index + Nx * Ny
				 alfa_N = alfa + 2 * Nx * Ny     ! Neighbour index in Y vector           
				 Del_PO = D_Poten( N_index  , B_index , 'o')				 
				 Del_PW = D_Poten( N_index  , B_index , 'w')	
				 				 
				 Tzao =	Transza(i,j,k,'o')
				 Tzaw =	Transza(i,j,k,'w')				 				 				 				 
	                          
				 ! d(fo(ijk))/d(po(ijk))
				 term = DTrza_DP(i,j,k,'o','c') * Del_PO + Tzao * ( - 1.0 - gama_der(B_index , 'o' ) * ( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa - 1) + 7
				 Jac(Jac_index) =Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa
				 				 							 
				 ! d(fw(ijk))/d(po(ijk))
				 term = DTrza_DP(i,j,k,'w','c') * Del_PW + Tzaw * ( - 1.0 - gama_der(B_index , 'w' ) * ( Zblock(N_index) - Zblock(B_index) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 7
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa	
				 
				 ! d(fo(ijk))/d(Sw(ijk))
				 term = DTrza_DS(i,j,k,'o','c') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa 
				 j_Jac(Jac_index) =	alfa + 1	
				 				 
				 ! d(fw(ijk))/d(Sw(ijk))
				 term = DTrza_DS(i,j,k,'w','c') * Del_PW + Tzaw * PCO_der(B_index) * ( 1.0 + gama_der( B_index , 'w' ) * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 8
				 Jac(Jac_index) = Jac(Jac_index) +  term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa + 1			
				 
				 ! d(fo(ijk))/d(po(ijk+1))  , P(ijk+1) = Y (alfa_N)
				 term = DTrza_DP(i,j,k,'o','a') * Del_PO + Tzao * ( 1.0 - gama_der( N_index , 'o' ) * ( Zblock( N_index ) - Zblock( B_index ) ) )  
				 Jac_index = 15 * ( alfa - 1) + 13
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N		
				 				 
				 ! d(fw(ijk))/d(po(ijk+1))
				 term = DTrza_DP(i,j,k,'w','a') * Del_PW + Tzaw * ( 1.0 - gama_der( N_index , 'w') * ( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 13
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N		
				 
				 ! d(fo(ijk))/d(Sw(ijk+1))  , Sw(ijk+1) = Y (alfa_N + 1)
				 term = DTrza_DS(i,j,k,'o','a') * Del_PO
				 Jac_index = 15 * ( alfa - 1) + 14
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa
				 j_Jac(Jac_index) =	alfa_N + 1						 
				 
				 ! d(fw(ijk))/d(Sw(ijk+1))
				 term = DTrza_DS(i,j,k,'w','a') * Del_PW + Tzaw * PCO_der(N_index) * ( - 1.0 + gama_der(N_index , 'w' ) *( Zblock( N_index ) - Zblock( B_index ) ) ) 
				 Jac_index = 15 * ( alfa + 1 - 1) + 14
				 Jac(Jac_index) = term
				 i_Jac(Jac_index) = alfa + 1
				 j_Jac(Jac_index) =	alfa_N + 1						 

			 EndIf
	         
			 ! accumulation terms :
			 X = C_Rock * ( pressure(B_index ) - P_ref )
			 poros_n  = phi(B_index) * (1 + X + 1 / 2.0 * X**2  ) 
	                  
			 Vb = Deltax(i) * Deltay(j) * Deltaz(B_index)
			 
			 ! d(fo(ijk))/d(Sw(ijk))
			 acc_term = Vb * poros_n / ( alfac * B_o(B_index) * d_time)         
			 Jac_index = 15 * ( alfa - 1) + 8
			 Jac(Jac_index) = Jac(Jac_index) +  acc_term
			 i_Jac(Jac_index) = alfa 
			 j_Jac(Jac_index) =	alfa + 1			 

			 ! d(fw(ijk))/d(Sw(ijk))
			 acc_term = - Vb * poros_n / ( alfac * B_w(B_index) * d_time) - Vb * poros_n * (Sw(B_index) ) / ( alfac * ( B_w(B_index) )**2 * d_time) * B_der( B_index , 'w' ) * PCO_der(B_index) 
			 Jac_index = 15 * ( alfa + 1 - 1) + 8
			 Jac(Jac_index) = Jac(Jac_index) +  acc_term
			 i_Jac(Jac_index) = alfa + 1
			 j_Jac(Jac_index) =	alfa + 1				 

			 ! d(fo(ijk))/d(Po(ijk))
			 acc_term = - Vb * phi(B_index) * ( C_Rock  + X * C_Rock ) * ( 1 - Sw(B_index) ) / ( alfac * B_o( B_index) * d_time)  + Vb * poros_n * ( 1 - Sw(B_index) ) / ( alfac * ( B_o(B_index) ) ** 2 * d_time) * B_der( B_index , 'o' )
			 Jac_index = 15 * ( alfa - 1) + 7
			 Jac(Jac_index) = Jac(Jac_index) +  acc_term
			 i_Jac(Jac_index) = alfa
			 j_Jac(Jac_index) =	alfa			 

			 ! d(fw(ijk))/d(Po(ijk))  , { d(fw)/d(po) = d(fw)/d(pw)}
			 acc_term = - Vb * phi(B_index) * ( C_Rock  + X * C_Rock ) * (Sw(B_index) ) / ( alfac * B_w(B_index) * d_time)  + Vb * poros_n * ( Sw(B_index) ) / ( alfac * ( B_w(B_index) ) ** 2 * d_time) * B_der( B_index , 'w' )
			 Jac_index = 15 * ( alfa + 1 - 1) + 7
			 Jac(Jac_index) = Jac(Jac_index) +  acc_term
			 i_Jac(Jac_index) = alfa + 1
			 j_Jac(Jac_index) =	alfa				 

		  EndDO
	   EndDO
	EndDO

	! Derivative of well equations

	Do i_well = 1 , Nwell
		W_index = Ntotal + i_well  ! Pwf position in pressure vector
		alfa_w = 2 * Ntotal + i_well 
		if ( W_type(i_well) == 'prod') then
		   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
			  WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2) - 1) * Nx +  ( i_k - 1) * Nx * Ny ! Well Block Index
			  B_alfa = 2 * ( WB_index ) - 1
	          
			  ! d(fo(ijk))/d(Po(ijk))
			  Well_term = - W_index_o(i_well,i_k) - W_index_o(i_well,i_k) * ( - B_der(WB_index,'o') / B_o(WB_index) - mu_der(WB_index,'o') / mu_f (pressure(WB_index),'o')) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k))
			  Jac_index = 15 * ( B_alfa - 1) + 7
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa
			  j_Jac(Jac_index) = B_alfa			  
	          
			  ! d(fo(ijk))/d(Pwf)          
			  Well_term = W_index_o(i_well,i_k) 
			  Jac_index = 15 * ( B_alfa - 1) + 15
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa
			  j_Jac(Jac_index) = alfa_w

			  ! d(fw(ijk))/d(Po(ijk)) = d(fw(ijk))/d(Pw(ijk))										
			  Well_term = - W_index_w(i_well,i_k) - W_index_w(i_well,i_k) * ( - B_der(WB_index,'w') / B_w(WB_index) ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k)) !! - mu_der(WB_index,'w') / mu_f ( W_press(WB_index) , 'w') 
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 7
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa + 1
			  j_Jac(Jac_index) = B_alfa				 

			  ! d(fw(ijk))/d(Pwf)          
			  Well_term = W_index_w(i_well , i_k) 
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 15
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa	+ 1
			  j_Jac(Jac_index) = alfa_w			  
	          
			  ! d(fo(ijk))/d(Sw(ijk))
			  Well_term = - Well_index_cons(i_well , i_k ) / ( mu_f (pressure(WB_index),'o') * B_f( pressure(WB_index),'o' ) ) * ( pressure( WB_index ) - pressure( W_index )  - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'o')
			  Jac_index = 15 * ( B_alfa - 1) + 8
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa 
			  j_Jac(Jac_index) = B_alfa + 1				  
	          
			  ! d(fw(ijk))/d(Sw(ijk))
			  Well_term = - Well_index_cons(i_well , i_k ) / ( mu_f ( pressure(WB_index),'w') * B_f( W_press(WB_index),'w' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'w') - W_index_w(i_well,i_k) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') ) * PCO_der( WB_index ) ! + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w')  !! Mu_w is constant
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 8
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa + 1
			  j_Jac(Jac_index) = B_alfa + 1					  
	         
			  if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! oil rate is given
				  ! d(fwell(ijk))/d(Po(ijk))
				  Well_term = W_index_o(i_well,i_k) + W_index_o(i_well,i_k) * ( - B_der(WB_index,'o') / B_f( pressure(WB_index),'o' ) - mu_der(WB_index,'o') / mu_f (pressure(WB_index),'o') ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				  Jac_index = 15 * ( alfa_w - 1) + 7
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa
	          				  
				  
				  ! d(fwell)/d(Sw(ijk))
				  Well_term =  Well_index_cons(i_well , i_k ) / ( mu_f (pressure(WB_index),'o') * B_f( pressure(WB_index),'o' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'o')
				  Jac_index = 15 * ( alfa_w - 1) + 8
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa + 1 				  

				  ! d(fwell)/d(Pwf)
				  Well_term = - W_index_o(i_well,i_k) 
				  Jac_index = 15 * ( alfa_w - 1) + 15
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = alfa_w
				  
			  ELSEif ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then ! BHP is given
				  ! d(fwell(ijk))/d(Po(ijk)) = 0
				  ! d(fwell(ijk))/d(Pwf) = 1
				  Jac_index = 15 * ( alfa_w - 1) + 15
				  Jac(Jac_index) = 1.000
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = alfa_w
			  ELSEif ( Schedule(i_well , 2 , sche_N) .EQ. 3 ) then ! Total Liquid Rate is Specified
				  ! d(fwell(ijk))/d(Po(ijk))
				  Well_term = W_index_o(i_well,i_k) + W_index_o(i_well,i_k) * ( - B_der(WB_index,'o') / B_f( pressure(WB_index),'o' ) - mu_der(WB_index,'o') / mu_f (pressure(WB_index),'o') ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) ) &
							  + W_index_w(i_well,i_k) + W_index_w(i_well,i_k) * ( - B_der(WB_index,'w') / B_f( W_press(WB_index),'w' ) ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				  Jac_index = 15 * ( alfa_w - 1) + 7
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa
	          				  
				  
				  ! d(fwell)/d(Sw(ijk))
				  Well_term =  Well_index_cons(i_well , i_k ) / ( mu_f (pressure(WB_index),'o') * B_f( pressure(WB_index),'o' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'o') &
      						  +  W_index_w(i_well,i_k) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * PCO_der( WB_index ) &
							  +  Well_index_cons(i_well , i_k ) / ( mu_f ( pressure(WB_index),'w') * B_f( W_press(WB_index),'w' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'w')							 
				  Jac_index = 15 * ( alfa_w - 1) + 8
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa + 1 				  

				  ! d(fwell)/d(Pwf)
				  Well_term = - W_index_o(i_well,i_k) - W_index_w(i_well,i_k)
				  Jac_index = 15 * ( alfa_w - 1) + 15
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = alfa_w
				  				  
			  EndIf
	          
		   EndDo
		Else   ! W_type(i_well) == 'inje'
		   Do i_K =  Well_k(i_well , 1) , Well_k(i_well ,2) 
			  WB_index = Well_loc(i_well , 1 ) + ( Well_loc(i_well , 2) - 1) * Nx +  ( i_k - 1) * Nx * Ny ! Well Block Index
			  B_alfa = 2 * ( WB_index ) - 1
	          
			  ! d(fo(ijk))/d(Po(ijk)) = 0, No oil is produced or injected 
	          

			  ! d(fw(ijk))/d(Po(ijk)) = d(fw(ijk))/d(Pw(ijk))
			  Well_term = - W_index_w(i_well,i_k) - W_index_w(i_well,i_k) * ( - B_der(WB_index,'w') / B_f( W_press(WB_index),'w' ) - mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) ) 			  
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 7
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa + 1
			  j_Jac(Jac_index) = B_alfa					  

	          
			  ! d(fo(ijk))/d(Sw(ijk)) = 0, No oil is produced or injected 

			  ! d(fw(ijk))/d(Sw(ijk))
			  !Well_term = - Well_index_cons(i_well , i_k ) / ( mu_f ( pressure(WB_index),'w') * B_f( W_press(WB_index),'w' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'w') - W_index_w(i_well,i_k) * (pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * PCO_der( WB_index ) 
!!			  For an Injection Well , k_rw is constant , so the derivative is Zero
			  Well_term = - W_index_w(i_well,i_k) * (pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * PCO_der( WB_index ) 
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 8
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa + 1
			  j_Jac(Jac_index) = B_alfa + 1					  
	          
			  ! d(fw(ijk))/d(Pwf)          
			  Well_term = W_index_w(i_well , i_k) 
			  Jac_index = 15 * ( B_alfa + 1 - 1) + 15
			  Jac(Jac_index) = Jac(Jac_index) +  well_term
			  i_Jac(Jac_index) = B_alfa	+ 1
			  j_Jac(Jac_index) = alfa_w			  
	          			  
	          
			  if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! injection rate is given
				  ! d(fwell(iwell))/d(Po(ijk))= d(fwell(iwell))/d(Pw(ijk))
				  Well_term = W_index_w(i_well,i_k) + W_index_w(i_well,i_k) * ( - B_der(WB_index,'w') / B_f( W_press(WB_index),'w' ) - mu_der(WB_index,'w') / mu_f ( W_press(WB_index),'w') ) * ( pressure(WB_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
				  Jac_index = 15 * ( alfa_w - 1) + 7
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa				  

				  ! d(fwell)/d(Sw(ijk))  
				  !Well_term = Well_index_cons(i_well , i_k ) / ( mu_f ( pressure(WB_index),'w') * B_f( W_press(WB_index),'w' ) ) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * Dkr_dS( WB_index , 'w') + W_index_w(i_well,i_k) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * PCO_der( WB_index )                            
				  Well_term =  W_index_w(i_well,i_k) * ( pressure( WB_index ) - pressure( W_index ) - Del_Pwf(i_well, i_k) ) * ( B_der(WB_index,'w') / B_f( W_press(WB_index),'w') + mu_der(WB_index,'w') / mu_f ( pressure(WB_index) , 'w') ) * PCO_der( WB_index )                            
				  Jac_index = 15 * ( alfa_w - 1) + 8
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = B_alfa + 1 				  
	                           
				  ! d(fwell(ijk))/d(Pwf)
				  Well_term = - W_index_w(i_well,i_k)                              
				  Jac_index = 15 * ( alfa_w - 1) + 15
				  Jac(Jac_index) = Jac(Jac_index) +  well_term
				  i_Jac(Jac_index) = alfa_w
				  j_Jac(Jac_index) = alfa_w			  				  
			  ELSEIF ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then ! BHP is given
				  ! d(fwell(ijk))/d(Po(ijk)) = 0
				  ! d(fwell(ijk))/d(Pwf) = 1
				  Jac_index = 15 * ( alfa_w - 1) + 15
				  Jac(Jac_index) = 1.000
			      i_Jac(Jac_index) = alfa_w
			      j_Jac(Jac_index) = alfa_w
			  EndIf	          
		   EndDo               
		Endif
	EndDo

	                         


	nz_num1 = 0
	Do ycounter = 1 , N_equ * 15
		if ( Jac(ycounter) /= 0 ) then
			nz_num1 = nz_num1 + 1
		EndIf
	Enddo

	Allocate ( ia1(nz_num1))
	Allocate ( ja1(nz_num1))
	Allocate ( a1(nz_num1))
!!	
!!	If (.NOT. Adjoint_Sol) Then 
	i = 1     
	Do ycounter = 1 , 15 * N_equ
		if ( Jac(ycounter) /= 0 ) then
			a1(i) = Jac(ycounter)
			ia1(i) = i_Jac(ycounter)
			ja1(i) = j_Jac(ycounter)
			i = i + 1
		EndIf
	EndDO
!!	Else
!!		i = 1     
!!		Do ycounter = 1 , 15 * N_equ
!!			if ( Jac(ycounter) /= 0 ) then
!!				a1(i) = Jac(ycounter)
!!				ia1(i) = i_Jac(ycounter)
!!				ja1(i) = j_Jac(ycounter)
!!				i = i + 1
!!			EndIf
!!		EndDO
!!	EndIf	    

	DeALLOCATE (Jac ,i_Jac ,j_Jac)	
	    

	EndSubroutine calc_Jac                                                                                                      
! ##########################################################################################################################################################
! ##########################################################################################################################################################

EndModule Mat_Calc