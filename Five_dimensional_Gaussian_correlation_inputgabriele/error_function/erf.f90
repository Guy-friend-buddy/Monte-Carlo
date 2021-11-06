	FUNCTION erf_s_sher(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: erf_s_sher
	erf_s_sher=gammp(0.5_sp,x**2)
	if (x < 0.0) erf_s_sher=-erf_s_sher
	END FUNCTION erf_s_sher


	FUNCTION erf_v_sher(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: erf_v_sher
	erf_v_sher=gammp(spread(0.5_sp,1,size(x)),x**2)
	where (x < 0.0) erf_v_sher=-erf_v_sher
	END FUNCTION erf_v_sher
