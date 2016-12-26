#ifndef ROM_ROBCODEC_H
#define ROM_ROBCODEC_H

#include "NodeDof6Buffer.h"

template <typename InFile, typename OutFile>
void transfer_rob(InFile &input, OutFile &output) {
  Rom::NodeDof6Buffer state_buffer(input.nodeCount());
  while (input.validCurrentState()) {
    const double header = input.currentStateHeaderValue();
    input.currentStateBuffer(state_buffer);
    output.stateAdd(state_buffer, header);
    input.currentStateIndexInc();
  }
}

template <typename InFile, typename OutFile>
void convert_rob(const std::string &inFilename, const std::string &outFilename) {
  InFile input(inFilename);
  OutFile output(outFilename, input.nodeIdBegin(), input.nodeIdEnd(), false);
  transfer_rob(input, output);
}

#endif /* ROM_ROBCODEC_H */
