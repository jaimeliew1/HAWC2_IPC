module IPC_types
implicit none    
integer, parameter :: mk = kind(1.0d0)

type IPC_var
  real*8, dimension(:,:) , allocatable ::  B, A, X , Y 
  real*8 :: KP
  integer :: N
end type IPC_Var

type(IPC_var), save :: IPC_array
   
contains
end module IPC_types
    