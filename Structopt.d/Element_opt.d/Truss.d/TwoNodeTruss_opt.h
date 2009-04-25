#ifndef _2nodetruss_OPT_HPP_
#define _2nodetruss_OPT_HPP_

#include <Element.d/Element.h>
#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

class TwoNodeTruss_opt : public TwoNodeTruss, public Element_opt 
{
public:
  TwoNodeTruss_opt(int* n) : Element(), TwoNodeTruss(n) {}
  Element *clone() { return new TwoNodeTruss_opt(*this); }
};

#endif
