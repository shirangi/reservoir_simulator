Subroutine con_chek

! this sub checks the converegence at each time step
!!
Use globalvar
Use Well
Use Trans
Use Mat_Calc
!-------------------------------------Variable Declaration Part-------------------------------------------------------
Implicit None

integer :: i_well , i, j , k , alfa, B_index 
real*8 :: Oil_division , W_division , Nomenator_new ,  Nomenator_past , Denomenator1 , Vb , MS , WNom_new, WNom_past , W_Denom
real*8 :: X	, WMS

!----------------
MS = 1d-5 ! convergence criteria
WMS = 1d-5
Nomenator_new = 0
Nomenator_past = 0
Denomenator1 = 0

WNom_new = 0
WNom_past = 0
W_Denom = 0 



Do k=1,Nz
   Do j=1,Ny
      Do i=1,Nx
         alfa = (k-1) * Nx * Ny + (j-1) * Nx + i
         Vb = Deltax(i)*Deltay(j)*Deltaz(alfa)
         X = C_Rock * ( pressure(alfa) - P_ref )
	     Nomenator_new = Nomenator_new + Vb * phi(alfa) * ( 1 +  X + 1.0 / 2.0 * X**2 ) * ( 1 - Sw(alfa) ) / ( alfac * B_o( alfa ) )          
         WNom_new = WNom_new + Vb * phi(alfa) * ( 1 +  X + 1.0 / 2.0 * X**2 ) * ( Sw(alfa) ) / ( alfac * B_w(alfa) )    
         X = C_Rock * ( past_p(alfa) - P_ref )
         Nomenator_past = Nomenator_past + Vb * phi(alfa) * ( 1 + X + 1.0 / 2.0 * X**2 )  * ( 1 - Past_Sw(alfa) ) / ( alfac * B_f( past_p(alfa) , 'o' ) )
         WNom_past = WNom_past + Vb * phi(alfa) * (1 + X + 1.0 / 2.0 * X**2 )  * ( Past_Sw(alfa) ) / ( alfac * B_f( past_pw(alfa), 'w' ) ) !  past_p(alfa) - PCO_Sw( Past_Sw(alfa))         
      EndDo
   EndDo
EndDo

call Calc_rate

Do i_well = 1, Nwell
    if ( W_type(i_well) == 'prod') then
        Denomenator1 = Denomenator1 + oilrate (i_well) * d_time
    EndIf
EndDo

Oil_division = abs(Nomenator_new  -  Nomenator_past) / Denomenator1



Do i_well = 1, Nwell
    W_Denom = W_Denom + waterrate (i_well) * d_time
EndDo

if (W_Denom .NE. 0.0) then
	W_division = abs(WNom_new  -  WNom_past) /abs(W_Denom)
Else 
	W_division = abs(WNom_new  -  WNom_past)
EndIF


call calc_RHS	 
!!    Open(1,DEFAULTFILE='.\sim_output\' , FILE='convergence.DAT', status='unknown', access = 'append')
			!Write(1,'(1X ," Time          BHP")')
			!write(1,'(I3,2x, F15.5)') time , pressure(Ntotal + iWell)
!!    write(1,'(F15.8,I3,<3>(2x,F15.8))') time, j_Ensemble, Oil_Division, W_division, maxFCN
!!    close(1)   
				write(*,*) '-------------------------------------------'
				write(*,*) ' Ensemble Number :' , j_Ens !!j_Ensemble				
				write(*,*) 't=', time ,',OMB=', Oil_Division , ',W MB=',W_Division , 'M-RHS=',maxFCN , ', iteration # ',i_iter
				write(*,*) '-------------------------------'

                

if (maxFCN < MRHS) then 
	If (Denomenator1 < 1.0) Then  !! There is no production in Res 
		if ( ( abs ( W_Division - 1 ) < OMB/100.0)) then 
			Converged = .True.
		EndIf
	Else	
		if ( ( abs ( Oil_Division - 1 ) < OMB)) then 
			Converged = .True.
		Elseif (maxFCN < MRHS) then 
			if ( ( abs ( Oil_Division - 1 ) < MOMB)) then 
				Converged = .True.
			EndIf
		EndIf
	EndIf
EndIf	

!if (maxFCN < 5.0d-1) then 
!	IF (d_time .LE. 1.0) then	
!		if ((i_iter > 10 )) then
!			Converged = .True.
!		EndIF		
!	EndIF
	if ((i_iter > M_ITER )) then
		Converged = .True.
	Endif
!ElseIF((i_iter > 30 )) then
!	Converged = .True.
!Endif

if (converged ) then
	write(*,*) '-------------------------------------------'       
	write(*,*) ' Ensemble Number :' , j_Ens !!! j_Ensemble 'MAP Estimate'
	write(*,*) '--C O N V E R G E D -- '
	write(*,*) '-------------------------------------------'       	   	
EndIf



EndSubroutine con_chek