# Makefile for limited memory bundle method
#FF = f77 +U77
#FF = f77
FF = g77 -Wall -g -fbounds-check

all:	tlmbm

tlmbm.o: tlmbm.f
	$(FF) -c tlmbm.f

lmbm.o: lmbm.f
	$(FF) -c lmbm.f

lmsub.o: lmsub.f
	$(FF) -c lmsub.f

matcal.o: matcal.f
	$(FF) -c matcal.f

tnsunc.o: tnsunc.f
	$(FF) -c tnsunc.f

tlmbm: tlmbm.o lmbm.o lmsub.o matcal.o tnsunc.o
	$(FF) -o tlmbm tlmbm.o lmbm.o lmsub.o matcal.o tnsunc.o

clean:	
	rm tlmbm tlmbm.o lmbm.o lmsub.o matcal.o tnsunc.o
