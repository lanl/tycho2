/*
    Assert.hh
    
    Macros for assert and insist.
    Assert can be turned on/off via ASSERT_ON.
    Insist is always checked.
*/

#ifndef __ASSERT_HH__
#define __ASSERT_HH__


namespace Assert
{
    void my_assert(const char *cond, const char *file, const int line);
    void my_insist(const char *cond, const char *msg, const char *file, const int line);
}



#ifndef ASSERT_ON
#define ASSERT_ON 1
#endif

#if ASSERT_ON
#define Assert(c) if (!(c)) Assert::my_assert( #c, __FILE__, __LINE__ );
#else
#define Assert(c) 
#endif

#define Insist(c,m) if (!(c)) Assert::my_insist( #c, m, __FILE__, __LINE__ );

#endif
