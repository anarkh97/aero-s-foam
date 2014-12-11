CXXFLAGS = -g -fPIC
AEROS_DIR = /home/avery/Codes/FEM

MatLaw.so: MatLaw.C
	g++ $(CXXFLAGS) -D_TEMPLATE_FIX_ -I$(AEROS_DIR) -c MatLaw.C ; g++ -shared MatLaw.o -o MatLaw.so

clean:
	rm -f MatLaw.o ; rm -f MatLaw.so ;

