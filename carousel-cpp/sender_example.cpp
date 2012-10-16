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

	kite::CarouselState cs;

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
