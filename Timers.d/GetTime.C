#include <sys/time.h>
#include <Timers.d/GetTime.h>

#ifdef WINDOWS
#include <sys/timeb.h>
#endif

#ifdef WINDOWS
double getTime()
{
	//return 0;
	timeb tb;
	ftime(&tb);
	double r = 1000.0*tb.time + tb.millitm/1000.0;
	return r;
}

#else
// function to get current time
double getTime() {
 timeval tp;
 struct timezone tz;
 gettimeofday(&tp, &tz);
 return 1000.0*tp.tv_sec + tp.tv_usec/1000.0;
}
#endif
