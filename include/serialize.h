#ifndef SERIALIZE_H
#define SERIALIZE_H
#include <stdint.h>
#include <memory.h>

static void* copy_to_double(void* src, void* dest)
{
  memcpy(dest, src, sizeof(double));
  return dest+sizeof(double); //FS this only works to compute the offset when copy to, not from...
}

static void* copy_to_double_array(void* src, void* dest, uint32_t size){
  memcpy(dest, src, size*sizeof(double));
  return dest+(size*sizeof(double));
}

static void* copy_to_int(void* src, void* dest)
{
  memcpy(dest, src, sizeof(int));
  return dest+sizeof(int);
}

static void* copy_to_long(void* src, void* dest)
{
  memcpy(dest, src, sizeof(long));
  return dest+sizeof(long);
}

static void* copy_from_double(void* src, void* dest)
{
  memcpy(dest, src, sizeof(double));
  return src+sizeof(double); //FS this only works to compute the offset when copy to, not from...
}

static void* copy_from_double_array(void* src, void* dest, uint32_t size){
  memcpy(dest, src, size*sizeof(double));
  return src+(size*sizeof(double));
}

static void* copy_from_int(void* src, void* dest)
{
  memcpy(dest, src, sizeof(int));
  return src+sizeof(int);
}

static void* copy_from_long(void* src, void* dest)
{
  memcpy(dest, src, sizeof(long));
  return src+sizeof(long);
}

#endif //SERIALIZE_H