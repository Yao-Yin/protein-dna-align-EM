#include "stdio.h"  
#include <omp.h>
  
int main(int argc,char * argv[])  
{  
    void* number =  0;  
    printf("%d\n",sizeof(&number));  
}  