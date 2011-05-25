#ifndef _ACCESS_H_
#define _ACCESS_H_

struct OffsetData {
  int first, last;
  double o[3];
};

class BoffsetAccessor{
 public:
  static int getNum(/*const*/ std::vector<OffsetData> & l, int i) 
    { 
      return ( (l[i]).last - (l[i]).first + 1);
    }
  static int getSize(const std::vector<OffsetData> & l)
    { return l.size(); }
  static int *getData(/*const*/ std::vector<OffsetData> & l, int i, int *nd)
    { 
      if(i<int(l.size()))
	{
	  int j = 0;
	  int * ndr;
	  if(nd != 0)
	    ndr = nd;
	  else
	    {
	      ndr = new int(getNum(l,i));
	      std ::cerr << "(W) warning : unefficient use of BoffsetAccessor" << std::endl;
	    }
          for(int k = (l[i]).first; k <= (l[i]).last; ++k,++j)
            ndr[j] = k;
	  return ndr;
	}
      else
	{
	  std::cerr << "BoffsetAccessor : offset out of range" << endl;
	}
      return(0);
    }
};

#endif
