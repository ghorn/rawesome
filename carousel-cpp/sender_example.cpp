// Quick and dirty cpp example that sends a kite protobuf message over zmq.

#include <string>
#include <iostream>
#include <fstream>
#include <stdint.h>

#include "kite.pb.h"

using namespace std;
int main(int argc, char **argv) 
{
	GOOGLE_PROTOBUF_VERIFY_VERSION;

	kite::Xyz xyz;
	xyz.set_x(2);
	xyz.set_y(2);
	xyz.set_z(2);

	kite::Dcm dcm;
	dcm.set_r11(1);
	dcm.set_r12(0);
	dcm.set_r13(0);
	dcm.set_r21(0);
	dcm.set_r22(1);
	dcm.set_r23(0);
	dcm.set_r31(0);
	dcm.set_r32(0);
	dcm.set_r33(1);

	kite::CarouselState cs;
	cs.mutable_kitexyz()->CopyFrom(xyz);
	cs.mutable_kitedcm()->CopyFrom(dcm);
	cs.set_delta(0.1);
	cs.set_rarm(1);
	cs.set_zt(-0.5);
	cs.set_w0(0);

	{
		// Write the kite state to disk.
		fstream output("sample_kite_message", ios::out | ios::trunc | ios::binary);
		if (!cs.SerializeToOstream(&output)) {
			cerr << "Failed to write sample_kite_message to disk." << endl;
			return -1;
		}
	}

	google::protobuf::ShutdownProtobufLibrary(); // optional
}
