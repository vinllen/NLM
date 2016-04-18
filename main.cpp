#include <map>

#include "process.h"
#include "nlm.h"
#include "unistd.h"

#define RENDER 0
using namespace std;

int main(int argc, char **argv)
{
    //int pic = atoi(argv[1]);
    //Process p(pic);
    
    int pic = -1;
    Process p(argv[1]); 

    p.readImage(); 
    pid_t pd = fork();
    if(pd == 0) {
        if (RENDER == 1)
            p.render();
    }
    else {
        puts("Input gaussian noise probability and Impluse probability:");
        double Gau = 0, Imp = 0, Rician = 0, level_p = 0;
        Gau = atoi(argv[2]);
        //Imp = atof(argv[3]);
        //Rician = atoi(argv[2]);
        p.addGaussianNoise(Gau);
        //p.addImpluseNoise(Imp);
        //p.addRicianNoise(Rician, level_p);
        //Gau = 15;
        p.writePic(p.getImageData(), "before_denoising.bmp");
        //Gau = level_p;
        puts("After add noise");

        pid_t pd2 = fork();
        if(pd2 == 0) {
            if (RENDER == 1)
                p.render();
        }
        else {
            NLM nlm(p.getImageData());
            int window, patch;
            window = 13, patch = 5;
            double s, sm, i, j, sm2;
            double m = 4 + 0.4 * Gau;
            if(Rician == 11)
                m = 22;
            else if(Rician == 9) 
                m = 20;
            else if(Rician == 7)
                m = 16;
            else if(Rician == 5)
                m = 14;
            else if(Rician == 3)
                m = 8;
            else if(Rician == 1)
                m = 4;
             
            string filename = "ans_classic";

            nlm.setWindow(window);
            nlm.setPatch(patch);
            nlm.setDelta(s, i, m, sm, j, sm2);
            puts("After set window and patch and delta");
            nlm.extendEdge();
            puts("After extend the edge");
            double final_fangcha = 0;
            //nlm.pretreatment(block_nr, final_fangcha);
            puts("After pretreatment");
            nlm.calcuW();
            puts("!!");
            if(UNBIAS)
                nlm.unbias(level_p);
            nlm.cutIllegal();
            p.setImageData(nlm.getImageData());
            p.caluSNR_PSNR(m, "ans.txt"); 
            //the parameter of jun,zhong,fangchar are just for print to file

            if (RENDER == 1)
                p.render();
        }
    }

    return 0;
}
