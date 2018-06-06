fast-com-mod:: main.o fast_mod_com_precompute_preinv.o geometric_interpolation.o
	gcc main.o fast_mod_com_precompute_preinv.o geometric_interpolation.o -o main -lflint -lmpfr -lgmp -lpthread
geometric_interpolation.o: geometric_interpolation.c geometric_interpolation.h
	gcc -c geometric_interpolation.c
fast_mod_com_precompute_preinv.o: fast_mod_com_precompute_preinv.c fast_mod_com_precompute_preinv.h geometric_interpolation.h
	gcc -c fast_mod_com_precompute_preinv.c
main.o: main.c fast_mod_com_precompute_preinv.h geometric_interpolation.h
	gcc -c main.c
clean:
	rm main.o
	rm geometric_interpolation.o
	rm fast_mod_com_precompute_preinv.o
