#include <stdio.h>
#define a(n) printf("%d,%d,%d,%d,",n+0, n+1, n+1, n+2)  // 2 bits
void b(unsigned n) { a(n+0),a(n+1),a(n+1),a(n+2); } // 4 bits
void c(unsigned n) { b(n+0),b(n+1),b(n+1),b(n+2); } // 6 bits
void d(unsigned n) { c(n+0),c(n+1),c(n+1),c(n+2); } // 8 bits
void e(unsigned n) { d(n+0),d(n+1),d(n+1),d(n+2); } //10 bits
void f(unsigned n) { e(n+0),e(n+1),e(n+1),e(n+2); } //12 bits
void g(unsigned n) { f(n+0),f(n+1),f(n+1),f(n+2); } //14 bits
void h(unsigned n) { g(n+0),g(n+1),g(n+1),g(n+2); } //16 bits
void i(unsigned n) { h(n+0),h(n+1),h(n+1),h(n+2); } //18 bits
void j(unsigned n) { i(n+0),i(n+1),i(n+1),i(n+2); } //20 bits
void k(unsigned n) { j(n+0),j(n+1),j(n+1),j(n+2); } //22 bits

int main()
{
    k(0);
}
