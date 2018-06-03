
Subroutine solve(N,AA,BB,XX)

use globalvar

implicit none

INTEGER :: N
Real*8, intent(out) :: XX(N)
Real*8, intent(in) :: BB(N) ,  AA(N,N)

REAL*8, ALLOCATABLE , DIMENSION(:,:):: L , U
REAL*8, ALLOCATABLE , DIMENSION(:):: bprim

integer :: i , j , k
real*8  :: sum


! to trasform an n*n matrix into LU

Allocate (L(N,N))
Allocate (U(N,N))
Allocate (bprim(N))


L = 0
U = 0
bprim = 0


Do i = 1 ,N
   L(i,1) = AA(i,1)
EndDo ! For i

!For j=1 To n
Do j = 1 ,N
 U(1,j) = AA(1,j)/ L(1,1)
EndDo ! for j


!pause

! For j = 2 To n
Do j = 2 , N
   !For i = j To n
   Do i = j , N
	sum = 0
	!For k = 1 To j-1
	Do k = 1 , j-1
	   Sum = sum + L(i,k)* U(k,j)
	EndDo ! for k
	L(i,j) = AA(i,j) - Sum
   EndDo ! for i

   U(j,j) = 1

   !For i = j + 1 To n
   Do  i = j + 1 , N
	sum = 0.0
	!For k = 1 To j-1
	Do k = 1 , j-1
	   Sum = Sum + L(j,k) * U (k,i)
	EndDo ! for k
	U(j,i) = (AA(j,i) - Sum) / L(j,j)
   EndDo !for i
EndDo ! for j





bprim(1) = BB (1) / L(1,1)

Do i = 2 , N
   Sum = 0 
   Do k = 1 , i-1
      sum = sum + L(i,k) * bprim ( k )
   End Do
   bprim(i) = (BB(i) - sum) / L(i,i)
EndDo

XX(N) = bprim (N)


Do j = N - 1 , 1 , -1
   sum = 0
   Do k = j+1 , N
      sum = sum + U(j,k) * XX(k)
   EndDo
   XX(j) = bprim(j) - sum
EndDo


DeAllocate (L)
DeAllocate (U)
DeAllocate (bprim)

End
