#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>

#define MISSING 99999
#define bipow(x,y) (((y)==0)? 1:(x))
#define logit(x) (log((x)/(1.0-x)))
#define inv_logit(x) ( 1.0 / ( 1.0 + exp(-(x)) ) )
#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))

#include "matrix.h"
#include "mathfunc.h"
#include "distribution.h"

long seed;

#include "datastruct.h" /*specifying data structure*/

/* Function getsize, below, determines the number of people in the
   population by counting the number of lines in the population input
   file. This function is used once at the start to assign the p_size
   variable. This will always be in increments of 2000 people*/


int get_size(FILE *in)
{
  int c, i=0;
  rewind(in);
  while((c=getc(in))!=EOF)
    if(c == '\n') i++;
  return i;
 }


// This program generates a pseudo-population of communities.
// There are five input parameters:
// (1) number of communities
// (2) number of people in each commuity
// (3) Last day of the epidemic
// (4) Is it a case-ascertained design, i.e., each community has a pre-fixed index case with illness onset on day 1
// (5) Is it a cluster randomization or individual randomization. 
// Only one time-independent is generated, and you can view it as treatment. Sometimes we also regard it as age group.
// We also genrate a time-dependent covariate from standard normal distribution.

int main(int argc, char *argv[])
{
  FILE *file;     
  PEOPLE *people; 
  PEOPLE *person;
  int h, i, j, k, l, n, m, r, t;
  int p_size, n_community, community_size;
  int day_epi_start, day_epi_stop;
  int case_ascertained, cluster_randomization;
  COMMUNITY *community;
  double prop_idx, prop_contact;
  void output(PEOPLE *people, int p_size, COMMUNITY *community, int n_community);

   
  /* This gets the random seed from the clock */

  seed = 123456789;
  mt_init(seed);

  if(argc != 6)
  {
    printf("Need 5 arguments to generate population\n");
    exit(0);
  }
  n_community = atoi(argv[1]);
  community_size = atoi(argv[2]);

  day_epi_start = 1;
  day_epi_stop = atoi(argv[3]);
  
  case_ascertained = atoi(argv[4]);
  cluster_randomization = atoi(argv[5]);
  /************************************************************************************************************/
  /************************************************************************************************************/
  /*** construct population structure ***/
  /************************************************************************************************************/
  /************************************************************************************************************/

  
  /* determine population file size */
  p_size = n_community * community_size;

  /* Allocate space for the people array. */
  people = (PEOPLE *)malloc((size_t) (p_size*sizeof(PEOPLE)));


  /* read the population file -- one line per person,
     so we loop over the size of the population */

  for(i=0,person=people; i < p_size; i++,person++)
  {
      person->id = i;
      person->community = i / community_size;
      make_1d_array_double(&person->time_ind_covariate, 1, 0.0);
      make_2d_array_double(&person->time_dep_covariate, 1, day_epi_stop - day_epi_start + 1, 0.0);
      person->pre_immune = 0;
      person->infection = 0;
      person->symptom = 0;
      person->day_ill = MISSING;
      person->exit = 0;
      person->day_exit = MISSING;
      person->idx = 0;
      person->u_mode = 0;
      person->q_mode = 0;
      person->ignore = 0;
      person->weight = 1;
  }
   
  /************************************************************************************************************/
  /************************************************************************************************************/
  /*** construct family structure ***/
  /************************************************************************************************************/
  /************************************************************************************************************/

  /* make space for the family numbers arrays */
  community = (COMMUNITY *)malloc((size_t) (n_community * sizeof(COMMUNITY)));

  /* Now record member id for each family. This is for the convenience
     of assigning quarantine to family members at the end of each day */

  for(h=0; h<n_community; h++) /* set size of each family to 0 */
  {
     community[h].size = 0;
     community[h].member = (int *)malloc((size_t) (community_size*sizeof(int)));
     community[h].day_epi_start = day_epi_start;
     community[h].day_epi_stop = day_epi_stop;
     community[h].day_last_followup = day_epi_stop;
     community[h].c2p_group = 0;
  }

  for(i=0,person=people;i<p_size;i++,person++) /* calculate family size */
  {
     h = person->community;
     community[h].member[community[h].size] = i;
     community[h].size++;
  }
  
  if(case_ascertained == 1) // set the first person in each community as th eindex case.
  {
     for(h=0; h<n_community; h++)
     {
        i = community[h].member[0];
        people[i].idx = 1;
        people[i].infection = 1;
        people[i].symptom = 1;
        people[i].day_ill = 1;
     }
  }

  // generate a time-independent covariate, which is treatment assignment in this  example.
  // Index cases and household contacts are randomized separately. For cluster randomization,
  // in the same household, all index cases receive the same treatment and all contacts receive the same treatment,
  // but the treatment can differ between index cases and contacts. 
  // For individual randomization, each individual is randomized separately. 
  prop_idx = 0.8;
  prop_contact = 0.5;
  if(cluster_randomization == 1)
  {
     for(h=0; h<n_community; h++)
     {
        k = (runiform(&seed) > prop_idx);
        l = (runiform(&seed) > prop_contact);
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].idx == 1) people[i].time_ind_covariate[0] = k;
           else people[i].time_ind_covariate[0] = l; 
        }   
     }
  }
  else
  {
     for(h=0; h<n_community; h++)
     {
        for(j=0; j<community[h].size; j++)
        {
           i = community[h].member[j];
           if(people[i].idx == 1) 
              people[i].time_ind_covariate[0] = (runiform(&seed) > prop_idx);
           else people[i].time_ind_covariate[0] = (runiform(&seed) > prop_contact); 
        }   
     }
  }
       
  // generate a time-dependent covariate using standard normal
  for(i=0; i < p_size; i++)
  {
      for(t=day_epi_start; t<=day_epi_stop; t++)
      {
         people[i].time_dep_covariate[0][t-1] = rnorm(0.0, 1.0, &seed);
      }
  }    

  output(people, p_size, community, n_community);

  for(i=0; i < p_size; i++)
  {
     free(people[i].time_ind_covariate);
     free_2d_array_double(people[i].time_dep_covariate);
  }
  free(people);
  for(h=0; h<n_community; h++)
  {
     free(community[h].member);
  }
  free(community);
}


void output(PEOPLE *people, int p_size, COMMUNITY *community, int n_community)
{
   int h, i, j, k, t;
   FILE *file;

   // generate input files for transtat
   file = fopen("pop.dat", "w");
   for(i=0; i < p_size; i++)
   {
      fprintf(file,"%8d  %8d  %1d  %1d  %1d  %8d  %1d  %8d  %1d  %2d  %2d  %f %1d\n",
              people[i].id, 
              people[i].community, 
              people[i].pre_immune, 
              people[i].infection, 
              people[i].symptom, 
              people[i].day_ill, 
              people[i].exit, 
              people[i].day_exit, 
              people[i].idx, 
              people[i].u_mode,
              people[i].q_mode,
              people[i].weight,
              people[i].ignore);
   }
   fclose(file);
   file = fopen("community.dat", "w");
   for(h=0; h < n_community; h++)
   {
      fprintf(file,"%8d  %8d  %8d  %8d  %8d\n",
                 h, community[h].day_epi_start, community[h].day_epi_stop, community[h].day_last_followup, community[h].c2p_group);
   }
   fclose(file);


   file = fopen("time_ind_covariate.dat", "w");
   for(i=0; i < p_size; i++)
   {
      fprintf(file,"%8d  %f\n", i, people[i].time_ind_covariate[0]);
   }
   fclose(file);

   file = fopen("time_dep_covariate.dat", "w");
   for(i=0; i < p_size; i++)
   {
      h = people[i].community;
      for(t=community[h].day_epi_start; t<=community[h].day_epi_stop; t++)
      fprintf(file,"%8d  %8d  %f\n", i, t, people[i].time_dep_covariate[0][t-1]);
   }
   fclose(file);
   /*
   file = fopen("c2p_contact.dat", "w");
   for(h=0; h<n_community; h++)
      fprintf(file,"%5d  %5d  %5d  %5d  %5.3f  %5d\n", h, community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0, 0);
   fclose(file);
   
   
   file = fopen("p2p_contact.dat", "w");
   for(h=0; h < n_community; h++)
      fprintf(file,"%5d  %5d  %5d  %5d  %5.3f  %5d\n", h, community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0, 0);
   fclose(file);
   */
}
