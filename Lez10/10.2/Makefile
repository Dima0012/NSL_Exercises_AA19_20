#LIBS:=`root-config --libs`
#INCS:=`root-config --cflags`

CC:= mpicxx
CFLAGS:= -Wall -O3 --std=c++11

HEADERS:= population.h individual.h random.h
OBJECTS:= main.o random.o individual.o population.o
MAIN:= main.exe
 
${MAIN}: ${OBJECTS}
	${CC} ${OBJECTS} -o $@ ${CFLAGS} #${INCS} ${LIBS}


%.o: %.cxx ${HEADERS}
	${CC} -c $< -o $@ ${CFLAGS} #${INCS} ${LIBS}

clean:
	rm *.o ${MAIN} seed.out

clean_data:
	rm Data/*.out
