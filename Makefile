OFILES = main.o ./source/ELASTOPLASTIC.o ./source/INTERPOLATION.o ./source/TOOLS.o ./source/PARTICLES.o ./source/GRID.o ./source/SIMULATION_PARAMETERS.o ./source/FILE_IO.o 
TARGET = mpm
CC = g++
CFLAGS = -Wall -pedantic -std=c++11 -O2

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@
	
clean:
	rm -f $(OFILES) $(TARGET)
	rm ./output/particle_position/*.txt
	rm ./output/particle_velocity/*.txt
	rm ./output/cauchy_stress/*.txt
	rm ./output/deformation_gradient/*.txt
	rm ./output/elastic_deformation_gradient/*.txt
	rm ./output/plastic_deformation_gradient/*.txt

# OUTPUT FROM "gcc -MM":

main.o: main.cpp include/ELASTOPLASTIC.h include/GRID.h \
  include/INTERPOLATION.h include/PARTICLES.h \
  include/SIMULATION_PARAMETERS.h include/TOOLS.h include/FILE_IO.h
ELASTOPLASTIC.o: source/ELASTOPLASTIC.cpp \
  source/../include/ELASTOPLASTIC.h
INTERPOLATION.o: source/INTERPOLATION.cpp \
  source/../include/INTERPOLATION.h
TOOLS.o: source/TOOLS.cpp source/../include/PARTICLES.h \
  source/../include/SIMULATION_PARAMETERS.h \
  source/../include/INTERPOLATION.h source/../include/TOOLS.h
FILE_IO.o: source/FILE_IO.cpp source/../include/FILE_IO.h
PARTICLES.o: source/PARTICLES.cpp source/../include/INTERPOLATION.h \
  source/../include/PARTICLES.h \
  source/../include/SIMULATION_PARAMETERS.h source/../include/TOOLS.h
GRID.o: source/GRID.cpp source/../include/GRID.h \
  source/../include/INTERPOLATION.h source/../include/ELASTOPLASTIC.h
SIMULATION_PARAMETERS.o: source/SIMULATION_PARAMETERS.cpp \
  source/../include/SIMULATION_PARAMETERS.h

