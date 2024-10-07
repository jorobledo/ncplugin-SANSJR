#include <iostream>
#include <vector>

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
    int n_points{1000};
    std::vector<double> q_vec(n_points);
    double sld{4.0};
    double solvent_sld{1.0};
    double radius{20.0};
    double length{400.0};
    double theta{60.0}, phi{60.0}, dtheta{0.0}, dphi{0.0};

    double step = (1.0 - 1e-3) / double(n_points - 1);
    double q=1e-3;

    QACRotation rotation;
    qac_rotation(&rotation, theta, phi, dtheta, dphi);
    std::cout << "#F1        F2        Iq" << std::endl;
    for (int i=0; i<n_points; ++i){
        q_vec[i] = q;
        // Anisotropic
        double qx = q * cos(theta) * cos(phi);
        double qy = q * cos(theta) * sin(phi);
        double qab=0.0, qc=0.0;
        qac_apply(&rotation, qx, qy, &qab, &qc);
        double Iq{Iqac(qab, qc, sld, solvent_sld, radius, length)};

        // Isotropic
        Fq(q, &F1, &F2, sld, solvent_sld, radius, length);
        F2 /= form_volume(radius, length);
        Iq /= form_volume(radius, length);
        std::cout <<  q << " " << F1 << " " << F2 << " " << Iq << std::endl;
        q+=step;
    }
    return 0;
}