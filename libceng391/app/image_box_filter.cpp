// Ertuğrul Demir
// 260201059
#include <cstdlib>
#include <iostream>

#include "ceng391/image.hpp"

using namespace std;
using ceng391::Image;

void exercise1_a(float sigma,std::string &filename);
void exercise1_b(float sigma,std::string &filename);
void exercise1_c(float sigma,std::string &filename);
void exercise2_a(std::string &filename);
void exercise2_b(std::string &filename);
void exercise3_b(float theta,std::string &filename);

int main(int argc, char** argv)
{
        if (argc < 2) {
                cerr << "Usage: " << argv[0] << " <filename> [<sigma>] [<theta>]" << endl;
                return EXIT_FAILURE;
        }

        float sigma = 3; // Default values
        int theta = 30;
        if (argc > 2)
                sigma = atoi(argv[2]);
        if (argc > 3)
                theta = atoi(argv[3]);

        exercise1_a (sigma, reinterpret_cast<string &>(argv[1]));
        exercise1_b (sigma, reinterpret_cast<string &>(argv[1]));
        exercise1_c (sigma, reinterpret_cast<string &>(argv[1]));
        exercise2_a (       reinterpret_cast<string &>(argv[1]));
        exercise2_b (       reinterpret_cast<string &>(argv[1]));
        exercise3_b (theta, reinterpret_cast<string &>(argv[1]));
        return EXIT_SUCCESS;
}

void exercise1_a(float sigma,std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->smooth_x(sigma);
    img->xsave_png("/tmp/smooth_x.png");
    delete img;
}

void exercise1_b(float sigma,std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->smooth_y(sigma);
    img->xsave_png("/tmp/smooth_y.png");
    delete img;
}

void exercise1_c(float sigma,std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->smooth(sigma);
    img->xsave_png("/tmp/smooth.png");
    delete img;
}

void exercise2_a(std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->deriv_x();
    img->xsave_png("/tmp/deriv_x.png");
    delete img;
}

void exercise2_b(std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->deriv_y();
    img->xsave_png("/tmp/deriv_y.png");
    delete img;
}

void exercise3_b(float theta,std::string &filename){
    Image *img =  Image::from_png(filename, false);
    img->rotate(theta,img);
    img->xsave_png("/tmp/rotate.png");
    delete img;
}

/// Local Variables:
/// mode: c++
/// compile-command: "make -C ../build image-box-filter"
/// End:
