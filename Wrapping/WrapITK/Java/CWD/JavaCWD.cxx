#include "JavaCWD.h"

#if defined(_WIN32) && (defined(_MSC_VER) || defined(__MINGW32__))
#include <direct.h>
#include <stdio.h>
#include <stdlib.h>

void JavaCWD::SetCWD(const char* dir)
{ 
  _chdir(dir);
}

const char* JavaCWD::GetCWD()
{ 
  static char buffer[4096];
  _getcwd(buffer, 4096);
  return buffer;
}
#else
#include <unistd.h>
void JavaCWD::SetCWD(const char* dir)
{ 
  chdir(dir);
}

const char* JavaCWD::GetCWD()
{ 
  static char buffer[4096];
  getcwd(buffer, 4096);
  return buffer;
}
#endif

#if !defined(_WIN32)
# include <dlfcn.h>

int JavaCWD::Load(const char* lib)
{
  return dlopen(lib, RTLD_GLOBAL|RTLD_NOW)? 1:0;
}
#else
int JavaCWD::Load(const char* lib)
{
  return 0;
}
#endif
