CXX=g++

TARGETS=sender_example
all: ${TARGETS} 

kite.pb.o:kite.pb.cc
	${CXX} kite.pb.cc -c -o kite.pb.o

sender_example:sender_example.o kite.pb.o
	${CXX} $< -o $@

clean:
	rm -f ${TARGETS} *.o 


