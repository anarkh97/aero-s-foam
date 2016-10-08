#ifndef _ATTRIB_H_
#define _ATTRIB_H_

// Element attribute structure
struct Attrib {
  int nele;
  int attr;
  int cmp_attr, cmp_frm;
  double cmp_theta;  // angle between element edge01 and material x axis
  Attrib() { nele = 0; attr = cmp_attr = cmp_frm = -1; cmp_theta = 0.0; }
};

#endif