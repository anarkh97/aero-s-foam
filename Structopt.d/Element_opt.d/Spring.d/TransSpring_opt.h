#ifndef _TRANSSPRLINK_OPT_HPP_
#define _TRANSSPRLINK_OPT_HPP_

#include <Element.d/Element.h>
#include <Element.d/Spring.d/TransSprlink.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

class TransSprlink_opt : public TransSprlink, public Element_opt 
{
public:
  TransSprlink_opt(int* n) : Element(), TransSprlink(n) {}
  Element *clone() { return new TransSprlink_opt(*this); }
};

#endif
