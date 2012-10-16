// Quick and dirty cpp example that sends a kite protobuf message over zmq.

#include <string>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <unistd.h>

#include <zmq.hpp>

#include "kite.pb.h"

using namespace std;
int main(int argc, char **argv) 
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	kite::Xyz xyz;
	kite::Dcm dcm;
	kite::CarouselState cs;

	xyz.set_x(2);
	xyz.set_y(2);
	xyz.set_z(2);
	dcm.set_r11(1);
	dcm.set_r12(0);
	dcm.set_r13(0);
	dcm.set_r21(0);
	dcm.set_r22(1);
	dcm.set_r23(0);
	dcm.set_r31(0);
	dcm.set_r32(0);
	dcm.set_r33(1);
	cs.mutable_kitexyz()->CopyFrom(xyz);
	cs.mutable_kitedcm()->CopyFrom(dcm);
	cs.set_delta(0.1);
	cs.set_rarm(1);
	cs.set_zt(-0.5);
	// To Add: Messages
	cs.set_w0(0);

#if 0
	{
		// Write the kite state to disk.
		fstream output("sample_kite_message", ios::out | ios::trunc | ios::binary);
		if (!cs.SerializeToOstream(&output)) {
			cerr << "Failed to write sample_kite_message to disk." << endl;
			return -1;
		}
	}
#else
	int i = 0;
	while(1)
	{
		cout << "Sending message of angle i = " << i << endl;
		// Send the kite state as a zmq message
		zmq::context_t context(1);
		zmq::socket_t socket(context,ZMQ_PUB);
		socket.bind("tcp://*:5563");
		string output;
		float delta_deg = (float)i / 360.0 * 2.0 * 3.1415;
		cs.set_delta(delta_deg);
		if (!cs.SerializeToString(&output)) {
			cerr << "Failed to serialize cs." << endl;
			return -1;
		}
		zmq::message_t message(output.size());
		memcpy((void *) message.data(), output.c_str(), output.size());
		socket.send(message);
		i += 10;
		if (i>350) i-=360;
		usleep(50*1000);
	}
#endif

	google::protobuf::ShutdownProtobufLibrary(); // optional
}
