#define PIPE
#ifndef PIPE

class Pipe {
public:
    Pipe();

    void setAngles(double, double, int); //min, max, num elements
    void setHeat(double, double, int);
    void setInletDiam(double, double, int);
    void setMassflow(double, double, int);

    void setGeometries(double*, double*, double*); // inlet diam, z, angle
    // calc surface area (2D), diameters (3D), cross areas (3D)

    ~Pipe(); 
private:
    // Input variables
    double *angle;
    double *heat_input;
    double *inlet_diameter;
    double *massflow;
    double *z;

    // Make these all double, triple pointers, etc so theyre easier to pass along to functions and use new keyword
    double **surfaceArea;
    double ***diameter;
    double ***crossArea;
    double *****quality;
};

#endif PIPE