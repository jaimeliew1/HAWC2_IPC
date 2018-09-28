!**************************************************************************************************
module IPC_HAWC2_main
use IPC_types
contains
subroutine init_ipc(array1,array2) bind(c, name='init_ipc')
   implicit none
   !DEC$ IF .NOT. DEFINED(__MAKEFILE__)
   !DEC$ ATTRIBUTES DLLEXPORT :: init_ipc
   !DEC$ END IF
 
  real*8, dimension(100) :: array1 , array2
  integer  N , l 
  real*8, dimension(:,:) , allocatable ::  B , A , X , Y
  real(mk) KP
 
  N =  array1(1)
  KP = array1(2)
  !allocate(IPC_array%X(3,N))
  !allocate(IPC_array%Y(3,N))
  allocate(IPC_array%X(3,N))
  allocate(IPC_array%Y(3,N))
  allocate(IPC_array%B(1,N))
  allocate(IPC_array%A(1,N))
 

  !allocate(B(1,N))
  !allocate(A(1,N))
  !allocate(X(1,N))
  !allocate(Y(1,N))
    IPC_array%X = 0
    IPC_array%Y = 0

 
   do l = 1,N
    !B(1,l) = array1(2+l)
    !A(1,l) = array1(N+2+l)
    IPC_array%B(1,l) = array1(2+l)
    IPC_array%A(1,l) = array1(N+2+l)
  enddo
 
 !IPC_array%X = X
 !IPC_array%Y = Y
 !IPC_array%B = B
 !IPC_array%A = A
 IPC_array%N = N
 IPC_array%KP = KP
 write(0,*) '! IPC Module with tracking. Jaime Liew 2018'
 


end subroutine init_ipc
!**************************************************************************************************
   subroutine update_ipc(array1, array2) bind(c, name='update_ipc')
   !DEC$ IF .NOT. DEFINED(__MAKEFILE__)
   !DEC$ ATTRIBUTES DLLEXPORT :: update_ipc
   !DEC$ END IF
   use IPC_types 
   implicit none
   real*8, dimension(100) :: array1 , array2
   real*8, dimension(:,:) , allocatable ::  X , Y , Xaux , Yaux
   real*8, dimension(1,3) :: theta , Xin , Theta_Out , Theta_aux
   integer :: N , N_var , q , w , Flag_debug
   real :: Xin_m , phase, amp, azim , deg2rad
  
   allocate(Xaux(3,N-1))
  
   N = IPC_array%N ;
   N_var = N-1 ;
 
   theta(1,1) = array1(1);
   theta(1,2) = array1(2);
   theta(1,3) = array1(3);
   Xin(1,1) = array1(4);
   Xin(1,2) = array1(5);
   Xin(1,3) = array1(6);
   azim = array1(7); 
   amp = array1(8);
   phase = array1(9); 

    Flag_debug = array1(10);
    
    Xin_m = (Xin(1,1) + Xin(1,2) + Xin(1,3)) / 3 
    Xin = Xin - Xin_m

    deg2rad = (3.14159265359/180) 
    Xin(1,1) = Xin(1,1) - amp*cos(azim*deg2rad + phase*deg2rad)
	Xin(1,2) = Xin(1,2) - amp*cos(azim*deg2rad - 120*deg2rad + phase*deg2rad)
	Xin(1,3) = Xin(1,3) - amp*cos(azim*deg2rad - 240*deg2rad  + phase*deg2rad)




   Xaux(1,:) = IPC_array%X(1,1:N-1)
  
   
IPC_array%X(1,2:N) = IPC_array%X(1,1:N_var) ;
IPC_array%X(2,2:N) = IPC_array%X(2,1:N_var) ;
IPC_array%X(3,2:N) = IPC_array%X(3,1:N_var) ;
IPC_array%X(:,1) = Xin(1,:) ;
 
IPC_array%Y(1,2:N) = IPC_array%Y(1,1:N_var) ;
IPC_array%Y(2,2:N) = IPC_array%Y(2,1:N_var) ;
IPC_array%Y(3,2:N) = IPC_array%Y(3,1:N_var) ;
IPC_array%Y(:,1) = 0 ;
   




! Calculation . Elementwise multiplication
 
if (Flag_debug.eq.1) then
write(0,*) '  Before Calculation '
write(0,*) 'IPC_array%B(1,w)', IPC_array%B
write(0,*) 'IPC_array%X', IPC_array%X
write(0,*) 'IPC_array%A', IPC_array%A
 
write(0,*) 'IPC_array%Y 1', IPC_array%Y(1,:)
write(0,*) 'IPC_array%Y 2', IPC_array%Y(2,:)
write(0,*) 'IPC_array%Y 3', IPC_array%Y(3,:)
endif
   do q = 1,3
       do w = 1,N

   IPC_array%Y(q,1) =  IPC_array%Y(q,1) + IPC_array%B(1,w)*IPC_array%X(q,w) - IPC_array%A(1,w)* IPC_array%Y(q,w)
  if (Flag_debug.eq.1)then
   write(0,*) ' IPC_array%Y(q,1)' ,  IPC_array%Y(q,1)
  
endif
       enddo
   enddo
! Output
Theta_aux(1,1) = IPC_array%Y(1,1)
Theta_aux(1,2) = IPC_array%Y(2,1)



Theta_aux(1,3) = IPC_array%Y(3,1)

 Theta_Out = theta - IPC_array%KP*Theta_aux
 
 
 
 

   
   
   
! Output
if (Flag_debug.eq.1)then
   write(0,*) '  After Calculation '
    write(0,*) '   '
      write(0,*) 'IPC_array%X(1,:)' , IPC_array%X(1,:)
      write(0,*) '   '
       write(0,*) 'IPC_array%Y(1,:)' , IPC_array%Y(1,:)
         write(0,*) '   '
        write(0,*) 'Theta_Out' , Theta_Out
         write(0,*) '   '
         write(0,*) 'IPC_array%KP', IPC_array%KP
         write(0,*) 'Theta_aux(1,1)',Theta_aux(1,1)
                 write(0,*) 'Theta_aux(1,2)',Theta_aux(1,2)
                         write(0,*) 'Theta_aux(1,3)',Theta_aux(1,3)
         write(0,*) 'Theta', theta
  pause
   endif
 

array2(1) = Theta_Out(1,1)
array2(2) = Theta_Out(1,2)
array2(3) = Theta_Out(1,3)



   end subroutine update_ipc
 
end module IPC_HAWC2_main
