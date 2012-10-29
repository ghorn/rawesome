// Quick and dirty cpp example that sends a kite protobuf message over zmq.

#include <string>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <unistd.h>

#include <zmq.hpp>
#include "zhelpers.hpp"

#include "kite.pb.h"

using namespace std;
int main(int argc, char **argv) 
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	kite::MultiCarousel mc;
	// Just pointers into mc
	kite::Xyz *xyz;
	kite::Dcm *dcm;
	kite::CarouselState *cs;

	mc.mutable_css();
	cs = mc.add_css();  // Combination rezize and return pointer to zeroth CarouselState (?)
	cs = mc.mutable_css(0);

	xyz = cs->mutable_kitexyz();
	dcm = cs->mutable_kitedcm();

	xyz->set_x(2);
	xyz->set_y(-0.5);
	xyz->set_z(0.5);
	dcm->set_r11(0); dcm->set_r12(1); dcm->set_r13(0);
	dcm->set_r21(-1); dcm->set_r22(0); dcm->set_r23(0);
	dcm->set_r31(0); dcm->set_r32(0); dcm->set_r33(1);
	cs->set_delta(0.1);
	cs->set_rarm(1);
	cs->set_zt(-0.05);
	cs->set_w0(0);

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
	zmq::context_t context(1);
	zmq::socket_t socket(context,ZMQ_PUB);
	socket.bind("tcp://*:5563");
	while(1)
	{
		cout << "Sending message of angle i = " << i << endl;
		// Send the kite state as a zmq message
		string output;
		float delta_deg = (float)i / 360.0 * 2.0 * 3.1415;
		cs->set_delta(delta_deg);
		//mc.add_css();
		//mc.mutable_css(0)->CopyFrom(cs); // runtime error
		if (!mc.SerializeToString(&output)) {
			cerr << "Failed to serialize mc." << endl;
			return -1;
		}
		s_sendmore(socket, "multi-carousel");
                s_send(socket, output);
		i += 10;
		if (i>350) i-=360;
		usleep(50*1000);
	}
#endif

	google::protobuf::ShutdownProtobufLibrary(); // optional
}
