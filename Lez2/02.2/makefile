#LIBS:=`root-config --libs`
#INCS:=`root-config --cflags`

CC:= g++
CFLAGS:= -Wall -O3 --std=c++11

HEADERS:= lib.h random.h walker.h
OBJECTS:= main.o lib.o random.o
MAIN:= main.exe
 
${MAIN}: ${OBJECTS}
	${CC} ${OBJECTS} -o $@ ${CFLAGS} #${INCS} ${LIBS}


%.o: %.cxx ${HEADERS}
	${CC} -c $< -o $@ ${CFLAGS} #${INCS} ${LIBS}

clean:
	rm *.o ${MAIN} seed.out Data/*.out
