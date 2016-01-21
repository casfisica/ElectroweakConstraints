#FLIBS=   /home/eduardo/lib32/libpacklib.a \
#         /home/eduardo/lib32/libmathlib.a 
#FLIBS=    -lpacklib -lmathlib
COMPILE.f=gfortran -c -C -g
LINK.f=gfortran -g
FFLAGS= -g

casobj  = main.f deltalib.f vmatrix.f #zprime.f


#          polylogs.o spence.o ffcli2.o ffxli2.o ffinit.o main.o zprime.o \
#          wwidth.o lep200.o taulife.o amu.o sumrules.o sin2th_plot.o


cas:    $(casobj)
	gfortran -o $@ $(FFLAGS) $(casobj) $(FLIBS)


clean:
	rm -rf *.o cas
cleanout:
	rm -rf *.out
