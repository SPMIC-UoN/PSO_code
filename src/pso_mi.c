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

#include "pso_mi.h"
#include <time.h>



//==============================================================
// calulate swarm size based on dimensionality
int pso_calc_swarm_size(int dim) {
  int size = 10. + 2. * sqrt(dim);
  return (size > PSO_MAX_SIZE ? PSO_MAX_SIZE : size);
}


//==============================================================
//          INERTIA WEIGHT UPDATE STRATEGIES
//==============================================================
// calculate linearly decreasing inertia weight
double calc_inertia_lin_dec_mi(int step, pso_settings_t_mi *settings_mi) {

  int dec_stage = 3 * settings_mi->steps / 4;
  if (step <= dec_stage)
    return settings_mi->w_min + (settings_mi->w_max - settings_mi->w_min) *  \
      (dec_stage - step) / dec_stage;
  else
    return settings_mi->w_min;
}



//==============================================================
//          NEIGHBORHOOD (COMM) MATRIX STRATEGIES
//==============================================================
// global neighborhood
void inform_global_mi(int *comm, double **pos_nb,
       double *pos_b, double *fit_b,
       double *gbest_mi, int improved,
       pso_settings_t_mi *settings_mi)
{
// printf("%d\n", 11111);
  int i;
  int j;
  // all particles have the same attractor (gbest)
  // copy the contents of gbest to pos_nb
  for (i=0; i<settings_mi->size; i++){

    for (j=0; j<settings_mi->dim; j++){
      pos_nb[i][j]=gbest_mi[j];
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

void inform_mi(int *comm, double *pos_nb, double *pos_b, double *fit_b,
      int improved, pso_settings_t_mi * settings_mi)
{
  int i, j;
  int b_n; // best neighbor in terms of fitness
  // for each particle
  for (j=0; j<settings_mi->size; j++) {
    b_n = j; // self is best
    // who is the best informer??
    for (i=0; i<settings_mi->size; i++)
      // the i^th particle informs the j^th particle
      if (comm[i*settings_mi->size + j] && fit_b[i] < fit_b[b_n])
        // found a better informer for j^th particle
        b_n = i;
    // copy pos_b of b_n^th particle to pos_nb[j]
    memmove((void *)&pos_nb[j*settings_mi->dim],
            (void *)&pos_b[b_n*settings_mi->dim],
            sizeof(double) * settings_mi->dim);

  }
}




// =============
// ring topology
// =============

// topology initialization :: this is a static (i.e. fixed) topology
void init_comm_ring_mi(int *comm, pso_settings_t_mi * settings_mi) {
  int i;
  // reset array
  memset((void *)comm, 0, sizeof(int)*settings_mi->size*settings_mi->size);

  // choose informers
  for (i=0; i<settings_mi->size; i++) {
    // set diagonal to 1
    comm[i*settings_mi->size+i] = 1;
    if (i==0) {
      // look right
      comm[i*settings_mi->size+i+1] = 1;
      // look left
      comm[(i+1)*settings_mi->size-1] = 1;
    } else if (i==settings_mi->size-1) {
      // look right
      comm[i*settings_mi->size] = 1;
      // look left
      comm[i*settings_mi->size+i-1] = 1;
    } else {
      // look right
      comm[i*settings_mi->size+i+1] = 1;
      // look left
      comm[i*settings_mi->size+i-1] = 1;
    }

  }

}




void inform_ring_mi(int *comm, double *pos_nb,
     double *pos_b, double *fit_b,
     double *gbest_mi, int improved,
     pso_settings_t_mi * settings_mi)
{

  // update pos_nb matrix
  inform_mi(comm, pos_nb, pos_b, fit_b, improved, settings_mi);

}

// ============================
// random neighborhood topology
// ============================
// void init_momm_random(int *comm, pso_settings_t * settings) {

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
//     init_momm_random(comm, settings);
//   inform_m(comm, pos_nb, pos_b, fit_b, improved, settings);

// }




//==============================================================
// return default pso settings
void pso_set_default_settings_mi(pso_settings_t_mi *settings_mi) {

  // set some default values
  settings_mi->dim = 4;//30;//22;
  // settings->x_lo = -20;
  // settings->x_hi = 20;
  settings_mi->goal = 1e-5;

  // settings->size = 8538;
  settings_mi->print_every = 10;
  settings_mi->steps = 500;
  settings_mi->c1 = 1.496;
  settings_mi->c2 = 1.496;
  settings_mi->w_max = 0.7298;
  settings_mi->w_min = 0.3;
  // settings->w_max = 0.000000001;
  // settings->w_min = 0.000000001;


  settings_mi->clamp_pos = 1;
  settings_mi->nhood_strategy = PSO_NHOOD_GLOBAL;
  settings_mi->nhood_size = 5;
  settings_mi->w_strategy = PSO_W_LIN_DEC;

  // settings->rng = NULL;
  // settings->seed = time(0);

}




//==============================================================
//                     PSO ALGORITHM
//==============================================================
void pso_solve_mi(pso_obj_fun_t_mi obj_fun_mi, void *obj_fun_params,
         pso_result_t_mi *solution_mi, pso_settings_t_mi *settings_mi, double *specs, double *B1inh, double *MT, double *APT, double *NOE)
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


  pos = malloc(settings_mi->size * sizeof(double*));
  vel = malloc(settings_mi->size * sizeof(double*));
  pos_b = malloc(settings_mi->size * sizeof(double*));
  pos_nb = malloc(settings_mi->size * sizeof(double*));


  int X;
  for(X=0;X<settings_mi->size;X++){
  pos[X] = (double *)malloc(settings_mi->dim * sizeof(double));
  vel[X] = (double *)malloc(settings_mi->dim * sizeof(double));
  pos_b[X] = (double *)malloc(settings_mi->dim * sizeof(double));
  pos_nb[X] = (double *)malloc(settings_mi->dim * sizeof(double));
  }

  fit = malloc(settings_mi->size * sizeof(double));
  fit_b = malloc(settings_mi->size * sizeof(double));


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
  comm = malloc(pow(settings_mi->size,2) * sizeof(int));
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
  switch (settings_mi->nhood_strategy)
    {
    case PSO_NHOOD_GLOBAL :
      // comm matrix not used
      inform_fun = inform_global_mi;
      break;
    case PSO_NHOOD_RING :
      init_comm_ring_mi((int *)comm, settings_mi);
      inform_fun = inform_ring_mi;
      break;
    // case PSO_NHOOD_RANDOM :
    //   init_momm_random((int *)comm, settings);
    //   inform_fun = inform_random;
    //   break;
    }

  // SELECT APPROPRIATE INERTIA WEIGHT UPDATE FUNCTION
  switch (settings_mi->w_strategy)
    {
      /* case PSO_W_mONST : */
      /*     calc_inertia_fun = calc_inertia_monst; */
      /*     break; */
    case PSO_W_LIN_DEC :
      calc_inertia_fun = calc_inertia_lin_dec_mi;
      break;
    }

  // INITIALIZE SOLUTION
  solution_mi->error_mi = DBL_MAX;

  // SWARM INITIALIZATION
  // for each particle
  for (i=0; i<settings_mi->size; i++) {
    // for each dimension
    for (d=0; d<settings_mi->dim; d++) {
      // generate two numbers within the specified range
      // a = settings->x_lo[d] + (settings->x_hi[d] - settings->x_lo[d]) * \
      //   gsl_rng_uniform(settings->rng);
      // b = settings->x_lo[d] + (settings->x_hi[d] - settings->x_lo[d]) *  \
      //   gsl_rng_uniform(settings->rng);
      srand(seed++);
      a = settings_mi->x_lo[d] + (settings_mi->x_hi[d] - settings_mi->x_lo[d]) * \
        rand() / 2147483647.0;
      srand(seed++);
      b = settings_mi->x_lo[d] + (settings_mi->x_hi[d] - settings_mi->x_lo[d]) *  \
        rand() / 2147483647.0;;
      // initialize position
      pos[i][d] = a;
      // best position is the same
      pos_b[i][d] = a;
      // initialize velocity
      vel[i][d] = (a-b) / 2.;
    }


    // update particle fitness
    fit[i] = obj_fun_mi(pos[i], settings_mi->dim, obj_fun_params, specs, B1inh, MT, APT, NOE);
    fit_b[i] = fit[i]; // this is also the personal best
    // update gbest??
    if (fit[i] < solution_mi->error_mi) {
      // update best fitness
      solution_mi->error_mi = fit[i];
      // copy particle pos to gbest vector
      int qq;
      for (qq=0; qq<settings_mi->dim; qq++){
        solution_mi->gbest_mi[qq]=pos[i][qq];
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
  for (step=0; step<settings_mi->steps; step++) {
    // update current step
    settings_mi->step = step;
    // update inertia weight
    // do not bother with calling a calc_w_monst function
    if (settings_mi->w_strategy)
      w = calc_inertia_fun(step, settings_mi);
    // check optimization goal
    if (solution_mi->error_mi <= settings_mi->goal) {
      // SOLVED!!
      if (settings_mi->print_every)
        printf("Goal achieved @ step %d (error=%.3e) :-)\n", step, solution_mi->error_mi);
      break;
    }

    // update pos_nb matrix (find best of neighborhood for all particles)
    inform_fun(comm, pos_nb, pos_b, fit_b, solution_mi->gbest_mi,
               improved, settings_mi);
    // the value of improved was just used; reset it
    improved = 0;
    // printf("%d\n", 70000);
    // update all particles
    for (i=0; i<settings_mi->size; i++) {
      // for each dimension

      for (d=0; d<settings_mi->dim; d++) {
        // calculate stochastic coefficients
        // rho1 = settings->c1 * gsl_rng_uniform(settings->rng);
        // rho2 = settings->c2 * gsl_rng_uniform(settings->rng);
        // printf("%d\n", 45678);
        srand(seed++);
        rho1 = settings_mi->c1 * rand() / 2147483647.0;
        srand(seed++);
        rho2 = settings_mi->c2 * rand() / 2147483647.0;
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
        if (settings_mi->clamp_pos) {
          if (pos[i][d] < settings_mi->x_lo[d]) {
            pos[i][d] = settings_mi->x_lo[d];
            vel[i][d] = 0;
          } else if (pos[i][d] > settings_mi->x_hi[d]) {
            pos[i][d] = settings_mi->x_hi[d];
            vel[i][d] = 0;
          }
        } else {
          // enforce periodic boundary conditions
          if (pos[i][d] < settings_mi->x_lo[d]) {

            pos[i][d] = settings_mi->x_hi[d] - fmod(settings_mi->x_lo[d] - pos[i][d],
                                              settings_mi->x_hi[d] - settings_mi->x_lo[d]);
            vel[i][d] = 0;

          } else if (pos[i][d] > settings_mi->x_hi[d]) {

            pos[i][d] = settings_mi->x_lo[d] + fmod(pos[i][d] - settings_mi->x_hi[d],
                                              settings_mi->x_hi[d] - settings_mi->x_lo[d]);
            vel[i][d] = 0;
          }
        }

      }

      // update particle fitness
      fit[i] = obj_fun_mi(pos[i], settings_mi->dim, obj_fun_params, specs, B1inh, MT, APT, NOE);
      // update personal best position?
      if (fit[i] < fit_b[i]) {
        fit_b[i] = fit[i];
        // copy contents of pos[i] to pos_b[i]
        // int ss;
        // for (ss=0; ss<settings->dim; ss++){
        //   pos[i][ss] = pos_b[i][ss];
        // }
        memmove((void *)&pos_b[i], (void *)&pos[i],
                sizeof(double) * settings_mi->dim);
      }

      // update gbest??
      if (fit[i] < solution_mi->error_mi) {
        improved = 1;
        // update best fitness
        solution_mi->error_mi = fit[i];
        // copy particle pos to gbest vector
        int qq1;
        for (qq1=0; qq1<settings_mi->dim; qq1++){
          solution_mi->gbest_mi[qq1]=pos[i][qq1];
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

    if (settings_mi->print_every && (step % settings_mi->print_every == 0))
      printf("Step %d (w=%.2f) :: min err=%.10f\n", step, w, solution_mi->error_mi);

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
