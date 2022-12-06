/* An implementation of the Particle Swarm Optimization algorithm

   Copyright 2010 Kyriakos Kentzoglanakis

   This program is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License version
   3 as published by the Free Software Foundation.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see
   <http://www.gnu.org/licenses/>.
*/



#include <time.h> // for time()
#include <math.h> // for cos(), pow(), sqrt() etc.
#include <float.h> // for DBL_MAX
#include <string.h> // for mem*
#include <stdlib.h>
#include <stdio.h>

// #include <gsl/gsl_rng.h>

#include "pso_di.h"
#include <time.h>



//==============================================================
// calulate swarm size based on dimensionality
int pso_calc_swarm_size_di(int dim) {
  int size = 10. + 2. * sqrt(dim);
  return (size > PSO_MAX_SIZE ? PSO_MAX_SIZE : size);
}


//==============================================================
//          INERTIA WEIGHT UPDATE STRATEGIES
//==============================================================
// calculate linearly decreasing inertia weight
double calc_inertia_lin_dec_di(int step, pso_settings_t_di *settings_di) {

  int dec_stage = 3 * settings_di->steps / 4;
  if (step <= dec_stage)
    return settings_di->w_min + (settings_di->w_max - settings_di->w_min) *  \
      (dec_stage - step) / dec_stage;
  else
    return settings_di->w_min;
}



//==============================================================
//          NEIGHBORHOOD (COMM) MATRIX STRATEGIES
//==============================================================
// global neighborhood
void inform_global_di(int *comm, double **pos_nb,
       double *pos_b, double *fit_b,
       double *gbest_di, int improved,
       pso_settings_t_di *settings_di)
{
// printf("%d\n", 11111);
  int i;
  int j;
  // all particles have the same attractor (gbest)
  // copy the contents of gbest to pos_nb
  for (i=0; i<settings_di->size; i++){

    for (j=0; j<settings_di->dim; j++){
      pos_nb[i][j]=gbest_di[j];
      // printf("%f ", pos_nb[i][j]);
    }
    // printf("\n");
    // memmove((void *)&pos_nb[i*settings->dim], (void *)&gbest,
    //         sizeof(double) * settings->dim);

    // int cc;
    // int dd;
    // printf("%d\n", i);
    // for (cc=0; cc<settings->dim; cc++){
    // printf("%.2f ", pos_nb[i*cc]);
    // }
    // printf("\n");
    // for (dd=0; dd<settings->dim; dd++){
    // printf("%.2f ", solution->gbest[dd]);
    // }
    // printf("\n");
  }
  // printf("%d\n", 99999);
}



// ===============================================================
// general inform function :: according to the connectivity
// matrix COMM, it copies the best position (from pos_b) of the
// informers of each particle to the pos_nb matrix

    // inform_fun(comm, pos_nb, pos_b, fit_b, solution->gbest,
    //            improved, settings);

void inform_di(int *comm, double *pos_nb, double *pos_b, double *fit_b,
      int improved, pso_settings_t_di * settings_di)
{
  int i, j;
  int b_n; // best neighbor in terms of fitness
  // for each particle
  for (j=0; j<settings_di->size; j++) {
    b_n = j; // self is best
    // who is the best informer??
    for (i=0; i<settings_di->size; i++)
      // the i^th particle informs the j^th particle
      if (comm[i*settings_di->size + j] && fit_b[i] < fit_b[b_n])
        // found a better informer for j^th particle
        b_n = i;
    // copy pos_b of b_n^th particle to pos_nb[j]
    memmove((void *)&pos_nb[j*settings_di->dim],
            (void *)&pos_b[b_n*settings_di->dim],
            sizeof(double) * settings_di->dim);

  }
}




// =============
// ring topology
// =============

// topology initialization :: this is a static (i.e. fixed) topology
void init_comm_ring_di(int *comm, pso_settings_t_di * settings_di) {
  int i;
  // reset array
  memset((void *)comm, 0, sizeof(int)*settings_di->size*settings_di->size);

  // choose informers
  for (i=0; i<settings_di->size; i++) {
    // set diagonal to 1
    comm[i*settings_di->size+i] = 1;
    if (i==0) {
      // look right
      comm[i*settings_di->size+i+1] = 1;
      // look left
      comm[(i+1)*settings_di->size-1] = 1;
    } else if (i==settings_di->size-1) {
      // look right
      comm[i*settings_di->size] = 1;
      // look left
      comm[i*settings_di->size+i-1] = 1;
    } else {
      // look right
      comm[i*settings_di->size+i+1] = 1;
      // look left
      comm[i*settings_di->size+i-1] = 1;
    }

  }

}




void inform_ring_di(int *comm, double *pos_nb,
     double *pos_b, double *fit_b,
     double *gbest_di, int improved,
     pso_settings_t_di * settings_di)
{

  // update pos_nb matrix
  inform_di(comm, pos_nb, pos_b, fit_b, improved, settings_di);

}

// ============================
// random neighborhood topology
// ============================
// void init_domm_random(int *comm, pso_settings_t * settings) {

//   int i, j, k;
//   // reset array
//   memset((void *)comm, 0, sizeof(int)*settings->size*settings->size);

//   // choose informers
//   for (i=0; i<settings->size; i++) {
//     // each particle informs itself
//     comm[i*settings->size + i] = 1;
//     // choose kappa (on average) informers for each particle
//     for (k=0; k<settings->nhood_size; k++) {
//       // generate a random index
//       j = gsl_rng_uniform_int(settings->rng, settings->size);
//       // particle i informs particle j
//       comm[i*settings->size + j] = 1;
//     }
//   }
// }
//
//
//
// void inform_random(int *comm, double *pos_nb,
//        double *pos_b, double *fit_b,
//        double *gbest, int improved,
//        pso_settings_t * settings)
// {


//   // regenerate connectivity??
//   if (!improved)
//     init_domm_random(comm, settings);
//   inform_d(comm, pos_nb, pos_b, fit_b, improved, settings);

// }




//==============================================================
// return default pso settings
void pso_set_default_settings_di(pso_settings_t_di *settings_di) {

  // set some default values
  settings_di->dim = 4;//30;//22;
  // settings->x_lo = -20;
  // settings->x_hi = 20;
  settings_di->goal = 1e-5;

  // settings->size = 8538;
  settings_di->print_every = 10;
  settings_di->steps = 500;
  settings_di->c1 = 1.496;
  settings_di->c2 = 1.496;
  settings_di->w_max = 0.7298;
  settings_di->w_min = 0.3;
  // settings->w_max = 0.000000001;
  // settings->w_min = 0.000000001;


  settings_di->clamp_pos = 1;
  settings_di->nhood_strategy = PSO_NHOOD_GLOBAL;
  settings_di->nhood_size = 5;
  settings_di->w_strategy = PSO_W_LIN_DEC;

  // settings->rng = NULL;
  // settings->seed = time(0);

}




//==============================================================
//                     PSO ALGORITHM
//==============================================================
void pso_solve_di(pso_obj_fun_t_di obj_fun_di, void *obj_fun_params,
         pso_result_t_di *solution_di, pso_settings_t_di *settings_di, double *specs, double *B1inh, double *MTi, double *APTi, double *NOEi, double *NOE2i)
{
  int free_rng = 0; // whether to free settings->rng when finished
  // Particles

  // MALLOCING
  double** pos;
  double** vel;
  double** pos_b;
  double* fit;
  double* fit_b;
  double** pos_nb;


  pos = malloc(settings_di->size * sizeof(double*));
  vel = malloc(settings_di->size * sizeof(double*));
  pos_b = malloc(settings_di->size * sizeof(double*));
  pos_nb = malloc(settings_di->size * sizeof(double*));


  int X;
  for(X=0;X<settings_di->size;X++){
  pos[X] = (double *)malloc(settings_di->dim * sizeof(double));
  vel[X] = (double *)malloc(settings_di->dim * sizeof(double));
  pos_b[X] = (double *)malloc(settings_di->dim * sizeof(double));
  pos_nb[X] = (double *)malloc(settings_di->dim * sizeof(double));
  }

  fit = malloc(settings_di->size * sizeof(double));
  fit_b = malloc(settings_di->size * sizeof(double));


  // NOT MALLOCING

  // double pos[settings->size][settings->dim]; // position matrix
  // double vel[settings->size][settings->dim]; // velocity matrix
  // double pos_b[settings->size][settings->dim]; // best position matrix
  // double fit[settings->size]; // particle fitness vector
  // double fit_b[settings->size]; // best fitness vector
  // pos[settings->size][settings->dim]; // position matrix
  // vel[settings->size][settings->dim]; // velocity matrix
  // pos_b[settings->size][settings->dim]; // best position matrix
  // fit[settings->size]; // particle fitness vector
  // fit_b[settings->size]; // best fitness vector

  // Swarm
  // double pos_nb[settings->size][settings->dim]; // what is the best informed
                                                // position for each particle
  int* comm;
  comm = malloc(pow(settings_di->size,2) * sizeof(int));
  //int comm[settings->size][settings->size]; // communications:who informs who
                                            // rows : those who inform
                                            // cols : those who are informed
  int improved; // whether solution->error was improved during
  // the last iteration

  int i, d, step;
  double a, b; // for matrix initialization
  double rho1, rho2; // random numbers (coefficients)
  double w; // current omega
  void (*inform_fun)(); // neighborhood update function
  double (*calc_inertia_fun)(); // inertia weight update function


  // CHECK RANDOM NUMBER GENERATOR
  // if (! settings->rng) {
  //   // initialize random number generator
  //   gsl_rng_env_setup();
  //   // allocate the RNG
  //   settings->rng = gsl_rng_alloc(gsl_rng_default);
  //   // seed the generator
  //   gsl_rng_set(settings->rng, settings->seed);
  //   // remember to free the RNG
  //   free_rng = 1;
  // }

  // initialise rng seed
  long seed = time(NULL);


  // SELECT APPROPRIATE NHOOD UPDATE FUNCTION
  switch (settings_di->nhood_strategy)
    {
    case PSO_NHOOD_GLOBAL :
      // comm matrix not used
      inform_fun = inform_global_di;
      break;
    case PSO_NHOOD_RING :
      init_comm_ring_di((int *)comm, settings_di);
      inform_fun = inform_ring_di;
      break;
    // case PSO_NHOOD_RANDOM :
    //   init_comm_random((int *)comm, settings);
    //   inform_fun = inform_random;
    //   break;
    }

  // SELECT APPROPRIATE INERTIA WEIGHT UPDATE FUNCTION
  switch (settings_di->w_strategy)
    {
      /* case PSO_W_CONST : */
      /*     calc_inertia_fun = calc_inertia_const; */
      /*     break; */
    case PSO_W_LIN_DEC :
      calc_inertia_fun = calc_inertia_lin_dec_di;
      break;
    }

  // INITIALIZE SOLUTION
  solution_di->error_di = DBL_MAX;

  // SWARM INITIALIZATION
  // for each particle
  for (i=0; i<settings_di->size; i++) {
    // for each dimension
    for (d=0; d<settings_di->dim; d++) {
      // generate two numbers within the specified range
      // a = settings->x_lo[d] + (settings->x_hi[d] - settings->x_lo[d]) * \
      //   gsl_rng_uniform(settings->rng);
      // b = settings->x_lo[d] + (settings->x_hi[d] - settings->x_lo[d]) *  \
      //   gsl_rng_uniform(settings->rng);
      srand(seed++);
      a = settings_di->x_lo[d] + (settings_di->x_hi[d] - settings_di->x_lo[d]) * \
        rand() / 2147483647.0;
      srand(seed++);
      b = settings_di->x_lo[d] + (settings_di->x_hi[d] - settings_di->x_lo[d]) *  \
        rand() / 2147483647.0;;
      // initialize position
      pos[i][d] = a;
      // best position is the same
      pos_b[i][d] = a;
      // initialize velocity
      vel[i][d] = (a-b) / 2.;
    }


    // update particle fitness
    fit[i] = obj_fun_di(pos[i], settings_di->dim, obj_fun_params, specs, B1inh, MTi, APTi, NOEi, NOE2i);
    fit_b[i] = fit[i]; // this is also the personal best
    // update gbest??
    if (fit[i] < solution_di->error_di) {
      // update best fitness
      solution_di->error_di = fit[i];
      // copy particle pos to gbest vector
      int qq;
      for (qq=0; qq<settings_di->dim; qq++){
        solution_di->gbest_di[qq]=pos[i][qq];
        // printf("%.1f " , solution->gbest[qq]);
      }
      // printf("\n%d\n", i);

    }
    // int p;
    // printf("%d\n", i);
    // for (p=0; p<30; p++){
    //   printf("%.1f ", pos[i][p]);
    // }
    //   printf("\n");
    // for (p=0; p<30; p++){
    //   printf("%.1f ", vel[i][p]);
    // }
    //   printf("\n");
    //   printf("%.1f ", fit[i]);
    //   printf("\n");
    //   printf("%.10f ", fit_b[i]);
    //   printf("\n");
    //   printf("%.10f ", solution->error);
    //   printf("\n");
    // for (p=0; p<30; p++){
    //   printf("%.1f " , solution->gbest[p]);
    // }
    //   printf("\n");
    //   printf("\n");
  }

  // initialize omega using standard value
  w = PSO_INERTIA;
  // RUN ALGORITHM
  for (step=0; step<settings_di->steps; step++) {
    // update current step
    settings_di->step = step;
    // update inertia weight
    // do not bother with calling a calc_w_dionst function
    if (settings_di->w_strategy)
      w = calc_inertia_fun(step, settings_di);
    // check optimization goal
    if (solution_di->error_di <= settings_di->goal) {
      // SOLVED!!
      if (settings_di->print_every)
        printf("Goal achieved @ step %d (error=%.3e) :-)\n", step, solution_di->error_di);
      break;
    }

    // update pos_nb matrix (find best of neighborhood for all particles)
    inform_fun(comm, pos_nb, pos_b, fit_b, solution_di->gbest_di,
               improved, settings_di);
    // the value of improved was just used; reset it
    improved = 0;
    // printf("%d\n", 70000);
    // update all particles
    for (i=0; i<settings_di->size; i++) {
      // for each dimension

      for (d=0; d<settings_di->dim; d++) {
        // calculate stochastic coefficients
        // rho1 = settings->c1 * gsl_rng_uniform(settings->rng);
        // rho2 = settings->c2 * gsl_rng_uniform(settings->rng);
        // printf("%d\n", 45678);
        srand(seed++);
        rho1 = settings_di->c1 * rand() / 2147483647.0;
        srand(seed++);
        rho2 = settings_di->c2 * rand() / 2147483647.0;
        // printf("%d\n", 87654);
        // update velocity
        // printf("%f\n", vel[i][d]);
        // printf("%f\n", pos[i][d]);
        // printf("%f\n", pos_b[i][d]);
        // printf("%ld\n", sizeof(pos_nb));
        // printf("%f\n", pos_nb[i][d]);
        vel[i][d] = w * vel[i][d] + \
          rho1 * (pos_b[i][d] - pos[i][d]) +  \
          rho2 * (pos_nb[i][d] - pos[i][d]);
          // printf("%d\n", 222);
        // update position
        pos[i][d] += vel[i][d];
        // clamp position within bounds?
        if (settings_di->clamp_pos) {
          if (pos[i][d] < settings_di->x_lo[d]) {
            pos[i][d] = settings_di->x_lo[d];
            vel[i][d] = 0;
          } else if (pos[i][d] > settings_di->x_hi[d]) {
            pos[i][d] = settings_di->x_hi[d];
            vel[i][d] = 0;
          }
        } else {
          // enforce periodic boundary conditions
          if (pos[i][d] < settings_di->x_lo[d]) {

            pos[i][d] = settings_di->x_hi[d] - fmod(settings_di->x_lo[d] - pos[i][d],
                                              settings_di->x_hi[d] - settings_di->x_lo[d]);
            vel[i][d] = 0;

          } else if (pos[i][d] > settings_di->x_hi[d]) {

            pos[i][d] = settings_di->x_lo[d] + fmod(pos[i][d] - settings_di->x_hi[d],
                                              settings_di->x_hi[d] - settings_di->x_lo[d]);
            vel[i][d] = 0;
          }
        }

      }

      // update particle fitness
      fit[i] = obj_fun_di(pos[i], settings_di->dim, obj_fun_params, specs, B1inh, MTi, APTi, NOEi, NOE2i);
      // update personal best position?
      if (fit[i] < fit_b[i]) {
        fit_b[i] = fit[i];
        // copy contents of pos[i] to pos_b[i]
        // int ss;
        // for (ss=0; ss<settings->dim; ss++){
        //   pos[i][ss] = pos_b[i][ss];
        // }
        memmove((void *)&pos_b[i], (void *)&pos[i],
                sizeof(double) * settings_di->dim);
      }

      // update gbest??
      if (fit[i] < solution_di->error_di) {
        improved = 1;
        // update best fitness
        solution_di->error_di = fit[i];
        // copy particle pos to gbest vector
        int qq1;
        for (qq1=0; qq1<settings_di->dim; qq1++){
          solution_di->gbest_di[qq1]=pos[i][qq1];
        }
        // memmove((void *)solution->gbest, (void *)pos[i],
        //         sizeof(double) * settings->dim);
      }
      // int p;
      // printf("\n\n%d\n", i);
      // for (p=0; p<settings->dim; p++){
      //   printf("%.1f ", pos[i][p]);
      // }
      //   printf("\n");
      // for (p=0; p<settings->dim; p++){
      //   printf("%.1f ", vel[i][p]);
      // }
      //   printf("\n");
      //   printf("%.1f ", fit[i]);
      //   printf("\n");
      //   printf("%.10f ", fit_b[i]);
      //   printf("\n");
      //   printf("%.10f ", solution->error);
      //   printf("\n");
      // for (p=0; p<settings->dim; p++){
      //   printf("%.1f " , solution->gbest[p]);
      // }
      //   printf("\n");
      //   printf("\n");
    }

    if (settings_di->print_every && (step % settings_di->print_every == 0))
      printf("Step %d (w=%.2f) :: min err=%.10f\n", step, w, solution_di->error_di);

  }
  // free(comm);
  // free(pos);
  // free(vel);
  // free(pos_b);

  // free(fit);
  // free(fit_b);
  // free(pos_nb);
  // free RNG??
  // if (free_rng)
    // gsl_rng_free(settings->rng);
  // free(rng);
}
