#linux:
CXXFLAGS = -O3 -fPIC

AEROS_DIR = /home/avery/Codes/FEM

#linux:
control.so: control.C
	g++ $(CXXFLAGS) -D_TEMPLATE_FIX_ -I$(AEROS_DIR) -c control.C ; g++ -shared control.o -o control.so

clean:
	rm -f control.o; rm -f control.so
