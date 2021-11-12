 
The objective of the code is to compare with the spectra calculation 

The input used here were provided by Gabriele. The input data comes from RANS + LES and no information on the scaling is avialable. 
Most possible all the length is scaled by chord and velocity by inlet velocity.  But cannot be confirmed. So present test case depends if or if not the scaling of spectra in gabriele's code is correct.  
The scaling used for the spectra of gabriele is 17961314271180.262  (= 2*np.pi*1./(velmean/1.43)*(1.225*16**2)**2/16*0.1356/4e-10)


Since the data provided was in matlab format I have saved it in .txt format
 file called ( input_gabriele), the data is saved in a double precision mode. The data is formated column wise i.e. 
each variable has a coloumn. 

Col 1 =  ypos/delta99  no dimension , this has been done since in the document gaussian.pdf equation 11 all the length scale has been normalised by delta (BL thickness)
Col 2 =  DU/Dy (mean shear) (1/s)
Col 3 = vertical velocity Std (m/s). 
Col 4 =    (1.4*4/sqrt(pi).*Lambda)/delta99----------- no dimension 

