/*#ifndef NUMTYPE_H
#define NUMTYPE_H
struct NumType {
NumType(double num);
NumType();
double fp;
int ep;
bool zf; // zeroflag
void rescaling();
void operator=(const NumType & other);
NumType operator*(const NumType & other);
NumType operator+(const NumType & other);
};
#endif*/