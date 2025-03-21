#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "constants.h"
#include <chrono>
#include "dynamics.h"
#include "osctomean.h"
#include "convert.h"
#include <eigen/Eigen/Dense>
#include "AstroCnvRef.h"
#include "calcnonlinearindx.h"
#include "uncertainityConvert.h"
 
#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace std;
using namespace matlab::data;

using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        double cst[15];

        cst[3] = inputs[0][0]; //Tmax
        cst[4] = inputs[0][1]; //m0
        cst[5] = inputs[0][2]; //mu
        cst[6] = inputs[0][3]; // J2
        cst[7] = inputs[0][4]; //Re

        const int method = inputs[0][5];

        cst[9] = inputs[0][6]; //params.LU;
        cst[10] = inputs[0][7]; //params.Area;
        cst[11] = inputs[0][8]; // params.Cd;
        cst[12] = inputs[0][9];//params.MU;
        cst[13] = inputs[0][13]; //param.drag

        const int Nnodes = inputs[0][10];
        const double tf = inputs[0][11];
        const double nu = inputs[0][12];
        const double coordSys = inputs[0][14];

        ArrayFactory f;

        // Cartesian 

        outputs[0] =f.createArray<double>({static_cast<unsigned long>(Nnodes),6});
        outputs[1] =f.createArray<double>({static_cast<unsigned long>(Nnodes),1});

        // MEE 

        outputs[2] =f.createArray<double>({static_cast<unsigned long>(Nnodes),6}); 
        outputs[3] =f.createArray<double>({static_cast<unsigned long>(Nnodes),1}); 

        // Geqoe

        outputs[4] =f.createArray<double>({static_cast<unsigned long>(Nnodes),6}); 
        outputs[5] =f.createArray<double>({static_cast<unsigned long>(Nnodes),1}); 

        // Ceqoe

        outputs[6] =f.createArray<double>({static_cast<unsigned long>(Nnodes),6}); 
        outputs[7] =f.createArray<double>({static_cast<unsigned long>(Nnodes),1}); 

        // Kepler 

        outputs[8] =f.createArray<double>({static_cast<unsigned long>(Nnodes),6}); 
        outputs[9] =f.createArray<double>({static_cast<unsigned long>(Nnodes),1}); 


        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 6);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), x0t(6), fx(6), fxbar(6);
  
        //DACE::AlgebraicMatrix<double> ;
        double statesGuess[6], xarray[6], constantPart[6];
        Eigen::MatrixXd Mat(6,6);
        AlgebraicVector<state_type> EigenV(6,0.0);
        AlgebraicVector<state_type> Geqoe_Eigen_DA(6);
        AlgebraicVector<state_type> Geqoesc(6);

        state_type nonlindx;


        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
        }

        
        int Jsum = 0;

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
        }

        // propagate x0, and output the result please.
        fx = RK78(6, x0, 0.0, tf, cst, DynamicsCart);

        if (k<1){
            fxbar = fx;
        }
        nonlindx = nonlinearindex(tf, fx,fxbar);

        for (unsigned int i = 0; i < 6; i++)
        {
            outputs[0][k][i] = cons(fx[i]);
        }

        outputs[1][k] = cons(nonlindx);


            //MEE OUTCOMES     
            //convert to mee
            x0t = eci2mee(x0, cst[5]);
            Geqoesc = convertU(tf, x0t);

            //propagate
            fx = RK78(6, mee2eci(Geqoesc, cst[5]), 0.0, tf, cst,
                      DynamicsCart);
            fx = eci2mee(fx, cst[5]);

            if (k<1){
                fxbar = fx;
            }
            nonlindx = nonlinearindex(tf, fx,fxbar);


            // get also the nonlinearity index.
            nonlindx = nonlinearindex(tf, fx,fxbar);

            for (unsigned int i = 0; i < 6; i++)
            { 
                outputs[2][k][i] = cons(fx[i]);
            }

             outputs[3][k]= cons(nonlindx);
       


            // Geqoe OUTCOMES 
            x0t     = RV2GEq(x0, cst);
            Geqoesc = convertU(tf, x0t);

            //propagate
            fx = RK78(6, GEq2RV(Geqoesc, cst), 0.0, tf, cst, 
                      DynamicsCart);
            fx = RV2GEq(fx, cst);

            if (k<1){
                fxbar = fx;
            }
            nonlindx = nonlinearindex(tf, fx,fxbar);


            for (unsigned int i = 0; i < 6; i++)
            { 
                outputs[4][k][i] = cons(fx[i]);
            }

            outputs[5][k] = cons(nonlindx);
     

            // Ceqoe OUTCOMES 
            x0t     = RV2CEq(x0, cst);
            Geqoesc = convertU(tf, x0t);

            //propagate
            fx = RK78(6, CEq2RV(Geqoesc, cst), 0.0, tf, cst, 
                      DynamicsCart);
            fx = RV2CEq(fx, cst);

               if (k<1){
                fxbar = fx;
            }
            nonlindx = nonlinearindex(tf, fx,fxbar);


            for (unsigned int i = 0; i < 6; i++)
            { 
                outputs[6][k][i] = cons(fx[i]);
                
            }

            outputs[7][k] = cons(nonlindx);

            // KEPLERIAN OUTCOMES
            x0t     = eci2Kep(x0,cst[5]);
            Geqoesc = convertU(tf, x0t);

            //propagate
            fx = RK78(6, Kep2eci(Geqoesc, cst[5]), 0.0, tf, cst, 
                      DynamicsCart);
            fx = eci2Kep(fx, cst[5]);

            if (k<1){
                fxbar = fx;
            }
            nonlindx = nonlinearindex(tf, fx,fxbar);


            for (unsigned int i = 0; i < 6; i++)
            { 
                outputs[8][k][i] = cons(fx[i]);
                
            }

            outputs[9][k] = cons(nonlindx);
    
        

        }
    

    };
