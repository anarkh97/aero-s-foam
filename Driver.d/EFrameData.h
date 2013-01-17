#ifndef _EFRAME_DATA_H_
#define _EFRAME_DATA_H_

// Element Frame Data structure
struct EFrameData {
    int elnum;
    EFrame frame;
    EFrameData *next;
};

// Node Frame Data structure
struct NFrameData {
    int elnum;
    EFrame frame;
    double origin[3];
    enum FrameType { Rectangular=0, Cylindrical, Spherical } type;
    NFrameData *next;
};

#endif
