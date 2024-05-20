#include <iostream>

extern "C" 
{
    #include "kernel_header.c"
    #include "lib/polevl.c"
    #include "lib/sas_J1.c"
    #include "lib/gauss150.c"
    #include "cylinder.c"
}

int main(int argc, char* argv[]){  

    double F1, F2;
    double q{0.1};
    double sld{1.0};
    double solvent_sld{2.0};
    double radius{20.0};
    double length{400.0};

    Fq(q, &F1, &F2, sld, solvent_sld, radius, length);
    std::cout << "Hello World!" << F1 << ", " << F2 << std::endl;
    return 0;
}