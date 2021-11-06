
 PROGRAM MAIN

 USE nrtype

 USE nr

 IMPLICIT NONE


 REAL(kind =8)  :: x,y

real(SP), external  :: erf_s_sher

x = 0.5d0
y = erf_s_sher(x) 

print *, "error function", y


end program main



