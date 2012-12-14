#linux:
CXXFLAGS = -O2 -fPIC
#sgi:
#CXXFLAGS = -g -64

#linux:
control.so: control.C
	g++ $(CXXFLAGS) -D_TEMPLATE_FIX_ -I/lustre/home/avery/hg/FEM -c control.C ; g++ -shared control.o -o control.so
#sgi:
#control.so: control.o
#	CC $(CXXFLAGS) -elf -shared -no_unresolved -rdata_shared -soname\
#           control.so \
#           -all control.o -o control.so -lm -lc
##	ld -elf -shared -no_unresolved -rdata_shared -soname control.so \
##           -all control.o -o control.so -lc

clean:
	rm -f control.o ; rm -f control.so
