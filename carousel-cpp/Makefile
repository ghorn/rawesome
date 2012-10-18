CXX=g++

TARGETS=sender_example
all: ${TARGETS} 

kite.pb.cc:../kite.proto
	protoc --cpp_out=. -I.. ../kite.proto

kite.pb.o:kite.pb.cc
	${CXX} kite.pb.cc -c -o kite.pb.o

sender_example.o:kite.pb.h

sender_example:sender_example.o kite.pb.o
	${CXX} $^ -o $@ -lprotobuf -lzmq

clean:
	rm -f ${TARGETS} *.o  *.pb.cc *.pb.h


