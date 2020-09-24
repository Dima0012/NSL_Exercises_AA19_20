#LIBS:=`root-config --libs`
#INCS:=`root-config --cflags`

CC:= g++
CFLAGS:= -Wall -O3 --std=c++11

HEADERS:= funzione.h  random.h integraleMC.h lib.h
OBJECTS:= main.o random.o integraleMC.o lib.o
MAIN:= main.exe
 
${MAIN}: ${OBJECTS}
	${CC} ${OBJECTS} -o $@ ${CFLAGS} #${INCS} ${LIBS}


%.o: %.cxx ${HEADERS}
	${CC} -c $< -o $@ ${CFLAGS} #${INCS} ${LIBS}

clean:
	rm *.o ${MAIN} seed.out Data/*.out
