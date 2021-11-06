
# fortran make file for P and L model using Gaussian correlation

OBJECTS = nrtype.o nr.o nrutil.o sobseq.o tridag.o locate.o spline.o splint.o moduleglobal.o gammln.o gser.o gcf.o gammp.o erf.o inversecdf.o brent_root.o contributiony.o contributionk.o main.o
.PHONY: clean

output.txt: main.exe
	./main.exe > output.txt

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe

%.o : %.f90
	gfortran -c $< 


wall_spectra_calculated.txt: main.exe   #
	./main.exe

plot: wall_spectra_calculated.txt       # 
	python compare_spectra.py


clean:
	rm -f $(OBJECTS) main.exe *.mod *.png wall_spectra_Panton.txt PDF.txt

