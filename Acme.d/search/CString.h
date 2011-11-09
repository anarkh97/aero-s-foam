// $Id$

#ifndef CString_h_
#define CString_h_
#include <cstddef>
#include <iostream>

#include "contact_assert.h"
#include "Contact_Defines.h"

const std::size_t C_NPOS = (std::size_t)-1;
// The maximum possible size for a string
const int C_STRING_DIGIT_LENGTH = 32;
// The maximum number of characters required to represent a number

/*****************************************************************************/
class CString 
/*****************************************************************************/
// This string class is intended to conform as closely as possible to
// the ANSI C++ standard library string class.
{
  public:
    CString(std::size_t n, char c);

    CString(const char *str = "");

    CString(const char *str, std::size_t n);

    CString(const CString &str, std::size_t pos = 0, std::size_t n = C_NPOS);

    ~CString();

    char  operator[](std::size_t pos) const;

    const char *data() const;
    const char *c_str() const;

    std::size_t size() const;

    std::size_t capacity() const;
    
    CString& append(const CString&);

/******************************* Implementation ******************************/

  private:
    struct my_Letter {
      std::size_t capacity;
      std::size_t size;
      std::size_t rc;
      char data[1];
    } *my_letter;

    void free(void);
};

inline std::size_t CString::size(void) const { return my_letter? my_letter->size : 0; }

inline std::size_t CString::capacity(void) const { 
  return my_letter? my_letter->capacity : 0; 
}

inline const char *CString::data(void) const 
{
  return my_letter->data;
}

inline const char *CString::c_str() const
{
  return my_letter->data;
}

inline char CString::operator[](std::size_t pos) const
{ 
  return my_letter->data[pos];
}

 
int operator==(const CString&, const char *);

std::ostream& operator<<(std::ostream &s, const CString &str);

CString operator<<(const CString &a, int b);

CString operator<<(const CString &a, Real b);

#endif // 
