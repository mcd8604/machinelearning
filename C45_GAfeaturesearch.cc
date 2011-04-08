#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include <stdio.h> // for P_tmpdir
#include <stdlib.h> // for mkstemp()
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
/* FOR SYSV Semaphores */
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/stat.h>
using namespace std;
int testcasenumber = 0;
const double tailwidth = 0.005;
int maxpoints = 1;
const int windowsize = 1;
int generations = 50;
int migration_interval = 10;
const int report_interval = 1;
int major_iteration = 0;
int popsize = 100;
const char verbose = 'y';
const char periodic_report = 'Y';
const int maxintweight = 20; //0,1,...,10
const int MAXTIES = 1000;
const int numbest = 1;
const int maxpopsize = 500;
const int MAXFEATURES = 50;
const int MAXTYPES = 500;
const int MAXBINS = 1024;
const int NUMCHILDREN = 2;
const int RATIO = 5; //  =  #cases/#bins
const int max_select = (maxpopsize*(maxpopsize+1))/2;
int select_vector[max_select];
int BINS;
int subsetsize;
int c_index, h_index;
int new_index;
int MINCASES = 5*RATIO;
double mutation_rate;
double crossover_rate = 0.95;
char const *buffer;
string linestring, fieldstring;
string test_case_label;
int numberoffields;
int numberselected;
int fieldend;
int numberoftypes = 0;
string typelist[MAXTYPES];
pid_t children_pids[NUMCHILDREN] = { 0 };
string children_popfiles[NUMCHILDREN];
string children_hashtables[NUMCHILDREN];
int typecount[MAXTYPES];
int testcases = 0;
double cost = 0.0;
int selected[MAXFEATURES];
double test_cases[windowsize][MAXFEATURES];
double examplebuffer[MAXFEATURES];
string class_label;
string class_labels[windowsize];
string predicted_label;
double field_value;
string popfilename;
string hashtablefilename;
int child_number=-1;

int semid = -1;

struct penalties
{
 string predictedname;
 double penalty;
 penalties * next;
};

struct classpenalty
{
 string classname;
 penalties * penaltylist;
 classpenalty * next;
};

struct fitnode
{
 int weights [MAXFEATURES];
 double fitness;
 fitnode * next;
};
const int buckets = 4723; //prime
fitnode* fhashtable[buckets] = { 0 };

struct newhashentry
{
  int bucket;
  fitnode const *entry;
  newhashentry *next;
};
newhashentry *newchildhashentries = NULL;


// The population is represented as an array of structs
struct chromosome
 {
  int ranking;
  double fitness;
  int weights[MAXFEATURES];
 };
chromosome bestN[maxpopsize];
chromosome newpop[maxpopsize];
chromosome population[maxpopsize];
chromosome overallbest;
chromosome previousbest;
chromosome children_populations[NUMCHILDREN][maxpopsize];
int rank_total;
double min_values[MAXFEATURES];
double max_values[MAXFEATURES];
double newmin[MAXFEATURES];
double newmax[MAXFEATURES];
double ranges[MAXFEATURES];
double totalcumulative = 0;

struct histostruct
{
 string type;
 double histogram [MAXFEATURES][MAXBINS];
 histostruct * next;
};
histostruct * aux;
histostruct * histolist;

struct testcase
{
 double featurearray[MAXFEATURES];
 string class_label;
 testcase * next;
};
 
string typenames[MAXTYPES];

double gettimeofday_dbl(void);
double absdiff(double n, double m);
int equal(int i, int j);
int select_chrom();
int hashfunction(int weights[]);
bool samestring(int a, int b);
int gauss(double mean, double stddev);
float gasdev();
void mate (int);
void rank(int, int);
void evaluate(fitnode * fhashtable[buckets], classpenalty * classpenaltylist,
               int member, testcase * testcaselist);
void evaluate_validation(classpenalty * classpenaltylist, chromosome overallbest,
              double validationpenalty[], int results, testcase * validationcaselist);
double penalty(classpenalty * classpenaltylist, string class_label, string
              predicted_label);
void updatecmatrix(string predicted_label, string class_label, 
                   int numberoftypes, string typelist[MAXTYPES],
                   int confmatrix[MAXTYPES][MAXTYPES]);

void updatehistos(double examplebuffer[],
             string  class_label,
             int numberoffields,
             double ranges[],
             double min_values[],
             double max_values[],
             histostruct * &histolist);

void updateglobalhistos(double examplebuffer[],
             int numberoffields,
             double ranges[],
             double min_values[],
             double max_values[],
             double globalhistos[MAXFEATURES][MAXBINS]);

string prediction(double test_cases[windowsize][MAXFEATURES],
                  histostruct * histolist,
                  int numberoffields,
                  int weights[MAXFEATURES],
                  double mins[MAXFEATURES],
                  double maxes[MAXFEATURES],
                  double ranges[MAXFEATURES]);

int read_popfile(char const * const popfilename,
                 chromosome * const population, 
                 int const popsize,
                 int const numselected);
int write_popfile(char const * const popfilename,
                  chromosome const * const population,
                  int const popsize,
                  int const numselected);
void append_to_hashtable(char const * const hashtablefilename,
                        int const numselected);
void write_new_hashtable_entries(char const * const hashtablefile,
                                 int const numselected);

int create_sema(void);
void destroy_sema(void);
void lock_sema(void);
void unlock_sema(void);

void generate_popfilenames (void);
void remove_popfiles (void);
void divvy_population(void);
bool spawn_children(void);
void setup_child(int child_num);
void wait_for_children(void);
void combine_children_results(double trainingpenalty[], bool newbestarray[], int results);

void generate_hashfilenames (void);
void remove_hashfiles (void);

int main(int argc, char* argv[])
{
ifstream examplefile(argv[1]);
ifstream testfile(argv[2]);
ifstream validationfile(argv[3]);
popfilename = "/dev/null";
srand(time(NULL));
//clock_t starttime, endtime;
double starttime, endtime;
//starttime = clock();
starttime = gettimeofday_dbl();
int i,j;
int results = 0;
double validationpenalty[1000];
double trainingpenalty[1000];
bool newbestarray[1000];
classpenalty * classpenaltylist = 0;
classpenalty * penaltyaux;
penalties * penaltyaux2;
cout<<"\nEnter class name. Enter # to stop.\n";
cin>>class_label;
while (class_label != "#")
 {
  penaltyaux = new classpenalty;
  penaltyaux -> classname = class_label;
  penaltyaux -> penaltylist = 0;
  cout<<"\nEnter labelled-as string. Enter # to stop."<<endl;
  cin>>predicted_label;
  while (predicted_label != "#")
   {
    penaltyaux2 = new penalties;
    penaltyaux2->predictedname = predicted_label;
    cout<<"\nEnter penalty.\n";
    cin>>penaltyaux2->penalty;
    penaltyaux2->next = (penaltyaux->penaltylist);
    (penaltyaux->penaltylist) = penaltyaux2;
    cout<<"\nEnter labelled-as string. Enter # to stop.\n";
    cin>>predicted_label;
   }
  penaltyaux -> next = classpenaltylist;
  classpenaltylist = penaltyaux;
  cout<<"\nEnter class name. Enter # to stop:\n";
  cin>>class_label;
 }
cout<<"\nPenalties entered:\n";
penaltyaux = classpenaltylist;
while (penaltyaux)
 {
  cout<<endl<<penaltyaux->classname<<": ";
  penaltyaux2 = penaltyaux->penaltylist;
  while (penaltyaux2)
   {
    cout<<penaltyaux2->predictedname<<", "<<penaltyaux2->penalty<<"; ";
    penaltyaux2 = penaltyaux2->next;
   }
  penaltyaux = penaltyaux->next;
 }
cout<<"\nEnter number of generations:\n";
cin>>generations;
cout<<"\nEnter population size:\n";
cin>>popsize;
cout<<"\nEnter migration interval:\n";
cin>>migration_interval;
cout<<"\nEnter number of fields:\n";
cin>>numberoffields;
cout<<"\nEnter number of fields in optimal subset:\n";
cin>>subsetsize;
if (subsetsize > numberoffields)
 subsetsize = numberoffields;
for (i=0; i<numberoffields; i++)
 {
  cout<<"\nEnter name of field "<<i<<": ";
  cin>>typenames[i];
 }
for (i=0; i<numberoffields; i++)
 cout<<typenames[i]<<endl;
cout<<"\nFor each field name listed, indicate with 'y' or 'n' "<<
 "whether or not to use the field.\n";
char response ='y';
for (i=0; i<numberoffields; i++)
 {
  cout<<"\n "<<typenames[i]<<" ?:\n";
  cin>>response;
  if (response == 'y')
   {
    numberselected++;
    selected[i] = true;
   }
  else
   selected[i] = false;
 }

//initialize hashtable
for (j=0; j<buckets; j++)
 fhashtable[j] = 0;

mutation_rate = 1.0/numberselected;

struct globalh
{
double globalhistos[MAXFEATURES][MAXBINS];
};

globalh * gptr = new globalh;

double freqsperfield[MAXFEATURES];
for (i=0; i<numberselected; i++)
 freqsperfield[i] = 0;
getline(examplefile, linestring);
  while (!examplefile.eof())
   {
    j = -1;
    fieldend = linestring.find_first_of(",",0);
    //erase record number
    linestring.erase(0,fieldend+1);
    fieldend = linestring.find_first_of(",",0);
    for (i=0; i<numberoffields; i++)
     {
      fieldend = linestring.find_first_of(",",0);
      fieldstring.assign(linestring, 0, fieldend);
      buffer = fieldstring.c_str();
      field_value = atof(buffer);
      if (selected[i])
       examplebuffer[++j] = field_value;
      linestring.erase(0,fieldend+1);
     }
   for (i=0; i<numberselected; i++)
    {
     if (freqsperfield[i] == 0)
      {
         freqsperfield[i] = 1;
         min_values[i] = examplebuffer[i];
         max_values[i] = examplebuffer[i];
       }
     else
      {
         if (examplebuffer[i] < min_values[i])
          min_values[i] = examplebuffer[i];
         if (examplebuffer[i] > max_values[i])
          max_values[i] = examplebuffer[i];
      }
    }
  getline(examplefile, linestring);
 }
for (i=0; i<numberselected; i++)
 ranges[i] = max_values[i] - min_values[i];
for (i=0; i<numberselected; i++)
 cout<<endl<<i<<": "<<max_values[i]<<" "<<min_values[i]
     <<" "<<ranges[i];


//create a global histogram
BINS = MAXBINS;
bool labelfound;
for (j=0; j<numberselected; j++)
 for (i=0; i<BINS; i++)
  gptr->globalhistos[j][i] = 0;
examplefile.clear();
examplefile.seekg(0,ios::beg);
getline(examplefile, linestring);
while (!examplefile.eof())
 {
  j = -1;
  fieldend = linestring.find_first_of(",",0);
  //erase record number
  linestring.erase(0,fieldend+1);
  for (i=0; i<numberoffields; i++)
   {
    fieldend = linestring.find_first_of(",",0);
    fieldstring.assign(linestring, 0, fieldend);
    buffer = fieldstring.c_str();
    field_value = atof(buffer);
    if (selected[i])
     examplebuffer[++j] = field_value;
    linestring.erase(0,fieldend+1);
   }
updateglobalhistos(examplebuffer, numberselected,
               ranges, min_values, max_values, gptr->globalhistos);
  getline(examplefile, linestring);
 }
int kk;
double total;
//convert to relative frequencies
  for (j=0; j<numberselected; j++)
   {
    total = 0.0;
    for (kk=0; kk<BINS; kk++)
     total += gptr->globalhistos[j][kk];
    for (kk=0; kk<BINS; kk++)
     gptr->globalhistos[j][kk] /= total;
   }
//for each feature, find tailwidth and 1-tailwidth values
 for (j=0; j<numberselected; j++)
  {
   i = 0;
   totalcumulative = 0;
   totalcumulative = gptr->globalhistos[j][i];
   while (totalcumulative <= tailwidth)
    totalcumulative += gptr->globalhistos[j][++i];
   //calculate newmin[j] from i,ranges,min_values,max_values
   newmin[j] = (i/(double)BINS) * ranges[j] + min_values[j];
cout<<"\nnewmin["<<j<<"] = "<<newmin[j];
   while (totalcumulative <= 1 - tailwidth)
    totalcumulative += gptr->globalhistos[j][++i];
   newmax[j] = (i/(double)BINS) * ranges[j] + min_values[j];
cout<<"\nnewmax["<<j<<"] = "<<newmax[j];
   min_values[j] = newmin[j];
   max_values[j] = newmax[j];
   ranges[j] = max_values[j] - min_values[j];
  }
delete gptr;

//create individual signal type histograms
histolist = 0;
examplefile.clear();
examplefile.seekg(0,ios::beg);
getline(examplefile, linestring);
while (!examplefile.eof())
 {
  j = -1;
  fieldend = linestring.find_first_of(",",0);
  //erase record number
  linestring.erase(0,fieldend+1);
  for (i=0; i<numberoffields; i++)
   {
    fieldend = linestring.find_first_of(",",0);
    fieldstring.assign(linestring, 0, fieldend);
    buffer = fieldstring.c_str();
    field_value = atof(buffer);
    if (selected[i])
     examplebuffer[++j] = field_value;
    linestring.erase(0,fieldend+1);
   }
  //get class name
  fieldend = linestring.find_first_of(",",0);
  fieldstring.assign(linestring, 0, fieldend);
  class_label = fieldstring;
  labelfound = false;
  for (int z=0; z<numberoftypes; z++)
   if (typelist[z] == class_label)
    {
     labelfound = true;
     typecount[z]++;
     break;
    }
  if (labelfound == false)
   {
    typelist[numberoftypes] = class_label;
    typecount[numberoftypes++] = 1;
   }
updatehistos(examplebuffer, class_label, numberselected,
               ranges, min_values, max_values, histolist);
  getline(examplefile, linestring);
 }
//convert to relative frequencies
int maxcount = typecount[0];
for (j=1; j<numberoftypes; j++)
 if (typecount[j] > maxcount)
  maxcount = typecount[j];
aux = histolist;
while (aux != 0)
 {
  for (j=0; j<numberselected; j++)
   for (kk=0; kk<BINS; kk++)
    aux -> histogram[j][kk] += 1.0/(maxcount+1);
  for (j=0; j<numberselected; j++)
   {
    total = 0.0;
    for (kk=0; kk<BINS; kk++)
     total += aux -> histogram[j][kk];
    for (kk=0; kk<BINS; kk++)
     aux -> histogram[j][kk] /= total;
   }
  aux = aux -> next;
 }
if (maxpoints > numberoftypes)
 maxpoints = numberoftypes;
  
//display histograms
//aux = histolist;
//while (aux != 0)
 //{
  //cout<<"\n\n\nType = "<<aux -> type<<endl;
  //for (j=0; j<numberoffields; j++)
   //{
    //cout<<"\nField number "<<j<<" :\n";
    //for (int k=0; k<BINS; k++)
     //cout<<aux -> histogram[j][k]<<"  ";
   //}
  //aux = aux -> next;
 //}

//create list of test cases
testcase * testcaselist = 0;
testcase * endlist = 0;
testcase * testptr;
getline(testfile, linestring);
while (!testfile.eof())
 {
  j= -1;
  testptr = new testcase;
  fieldend = linestring.find_first_of(",",0);
  //erase record number
  linestring.erase(0,fieldend+1);
  for (i=0; i<numberoffields; i++)
   {
    fieldend = linestring.find_first_of(",",0);
    fieldstring.assign(linestring, 0, fieldend);
    buffer = fieldstring.c_str();
    field_value = atof(buffer);
    if (selected[i])
     testptr->featurearray[++j] = field_value;
    linestring.erase(0,fieldend+1);
   }
  //get class name
  fieldend = linestring.find_first_of(",",0);
  fieldstring.assign(linestring, 0, fieldend);
  testptr->class_label = fieldstring;
  if (endlist)
   endlist->next = testptr;
  else
   testcaselist = testptr;
  testptr->next = 0;
  endlist = testptr;
  getline(testfile, linestring);
 }
//display testcaselist for debugging purposes
//cout<<endl<<"testcaselist: \n";
//testptr = testcaselist;
//while (testptr)
//{
 //for (j=0; j<numberselected; j++)
  //cout<<testptr->featurearray[j]<<" ";
 //cout<<endl;
 //cout<<testptr->class_label<<endl;
 //testptr = testptr->next;
//}

//create list of validation cases
testcase * validationcaselist = 0;
endlist = 0;
getline(validationfile, linestring);
while (!validationfile.eof())
 {
  j = -1;
  testptr = new testcase;
  fieldend = linestring.find_first_of(",",0);
  //erase record number
  linestring.erase(0,fieldend+1);
  for (i=0; i<numberoffields; i++)
   {
    fieldend = linestring.find_first_of(",",0);
    fieldstring.assign(linestring, 0, fieldend);
    buffer = fieldstring.c_str();
    field_value = atof(buffer);
    if (selected[i])
     testptr->featurearray[++j] = field_value;
    linestring.erase(0,fieldend+1);
   }
  //get class name
  fieldend = linestring.find_first_of(",",0);
  fieldstring.assign(linestring, 0, fieldend);
  testptr->class_label = fieldstring;
  if (endlist)
   endlist->next = testptr;
  else
   validationcaselist = testptr;
  testptr->next = 0;
  endlist = testptr;
  getline(validationfile, linestring);
 }
//display validationcaselist for debugging purposes
//cout<<endl;
//testptr = validationcaselist;
//while (testptr)
//{
 //for (j=0; j<numberselected; j++)
  //cout<<testptr->featurearray[j]<<" ";
 //cout<<endl;
 //testptr = testptr->next;
//}

//initialize population
int location, zzz;
for (i=0; i<popsize; i++)
 {
  for (j=0; j<numberselected; j++)
   population[i].weights[j] = 0;
  for (zzz=0; zzz<subsetsize; zzz++)
   {
    location = (rand()/100)%numberselected;
    while (population[i].weights[location] == 1)
     location = (rand()/100)%numberselected;
   population[i].weights[location] = 1;
  }
 }
//write population to file, comma separated weights
//write_popfile(argv[3], population, popsize, numberselected);

if (verbose == 'y')
{
 for (i=0; i<popsize; i++)
{
  cout<<endl;
  for (j=0; j<numberselected; j++)
   cout<<population[i].weights[j]<<"_";
 }
cout<<endl;
}

generate_popfilenames();
generate_hashfilenames();
create_sema();

overallbest.fitness = -50000000;
int gen;
int count;
int weight_value;

cout << "Generations: " << generations << endl;
cout << "Migration Interval: " << migration_interval << endl;

for (gen=0; gen<generations; gen += migration_interval)
{
 ++major_iteration;
 divvy_population();
 bool is_child_process = spawn_children();

 if (is_child_process)
 {
   for (int g = 0; g < migration_interval; ++g)
     {
       //evaluate current population
       for (i=0; i<popsize; i++)
         evaluate(fhashtable, classpenaltylist, i, testcaselist);
 
       //insert numbest best from previous population
       if (g > 1)
         for (i=0; i<numbest; i++)
           population[i] = bestN[i];

       //assign ranks, also fills out bestN with information
       // to use on the next iteration.
       rank(popsize, numbest); //use rank-based selection


       if ((verbose == 'y') && ((gen+g) % report_interval == 0))
         {
           lock_sema();

           // display ranked population, with fitness values
           cout<<endl<<"Generation "<<gen
               <<'+'<<g<<" child #"<<child_number
               <<" (pid: "<<getpid()<<"):";
           for (i=0; i<popsize; ++i)
             {
               cout<<endl;
               for (j=0; j<numberselected; ++j)
                 cout<<population[i].weights[j]<<'_';
               cout<<", "<<population[i].fitness;
             }
           cout<<endl;

           unlock_sema();
         }

       if (periodic_report == 'Y')
         if ((gen+g) % report_interval == 0)
           if (g+1 != migration_interval)
             {
               lock_sema();

               cout << endl <<
                 "generation "<<gen
                 <<'+'<<g<<" child #"<<child_number
                 <<" (pid: "<<getpid()<<") best so far: ";
               //cout << overallbest.fitness << ": ";
               cout << population[popsize-1].fitness << ": ";
               for (j=0; j<numberselected; ++j)
                 //cout << overallbest.weights[j] << " ";
                 cout << population[popsize-1].weights[j] << " ";
               cout<<endl;

               unlock_sema();
             }


       // If this is not the last generation, mate.
       if (g + 1 < migration_interval)
          {
           //form selection vector
           c_index=1;
           select_vector[0] = 0;
           for (i=1;i<popsize;i++)
            {
             if (samestring(i,i-1)) continue;
             count = population[i].ranking;
             for (j=1;j<=count;j++)
              {select_vector[c_index] = i; c_index++;}
            }
           //create new population
           new_index = -1;
           for (i=1; i<=popsize/2; i++) //assumes even population size
             mate(g);
              //selects from population[], populates newpop[]
           for (i=0; i<popsize; i++) //copy newpop[] into population[]
            population[i] = newpop[i];
         }

     }

    cerr << "Pid " << getpid() << " writing out " << popsize
         << " entries to " << popfilename << endl;
    write_popfile(popfilename.c_str(), population, popsize, numberselected);

    write_new_hashtable_entries(hashtablefilename.c_str(), numberselected);
    // Children don't make it out of this block of code.
    exit(0);
 }
 else
 {
   wait_for_children();
   combine_children_results(trainingpenalty, newbestarray, results);
   evaluate_validation(classpenaltylist, overallbest, validationpenalty, results, validationcaselist);
   results++;
 }

/* This output is useless, as things have not been ranked()
 * at this "pangea" level.
if ((verbose == 'y') && (gen % report_interval == 0))
{
  //display ranked population, with fitness values
  for (i=0; i<popsize; i++)
   {
    cout<<endl;
    for (j=0; j<numberselected; j++)
     cout<<population[i].weights[j]<<'_';
    cout<<", "<<population[i].fitness;
   }
  cout<<endl;
}

 if (periodic_report == 'Y')
   if (gen % report_interval == 0)
    if (gen != generations)
     {
      cout << endl <<
      "generation " << gen << " best so far: ";
      cout << overallbest.fitness << ": ";
      for (j=0; j<numberselected; j++)
       cout << overallbest.weights[j]<<" ";
      cout<<endl;
     }
*/
 if (results > 1)
  if (-1 * validationpenalty[results-1] > -1.1 * validationpenalty[results-2])
   {
    overallbest = previousbest;
    lock_sema();
    cout<<"\nNew validation penalty = "<<validationpenalty[results-1];
    cout<<"\nOld validation penalty = "<<validationpenalty[results-2];
    unlock_sema();
    break;
   }
    
}

rank(popsize, numbest);  // Rank the final results.

cout<<"\nFinal generation:\n";
 for (i=0; i<popsize; i++)
   {
    cout<<endl;
    for (j=0; j<numberselected; j++)
     cout<<population[i].weights[j]<<'_';
    cout<<", "<<population[i].fitness;
   }
  cout<<endl;

//cout<<"\nFinal generation fitness = "<<overallbest.fitness<<
cout<<"\nOverallbest fitness = "<<overallbest.fitness<<
    "\nweights = \n";
j=-1;
for (i=0; i<numberoffields; i++)
 if (selected[i])
  //cout<<typenames[i]<<": "<<overallbest.weights[++j]/(double)maxintweight<<endl;
  cout<<typenames[i]<<": "<<overallbest.weights[++j]<<endl;
j=-1;
for (i=0; i<numberoffields; i++)
 if (selected[i] && overallbest.weights[++j])
  cout<<"y\n";
 else
  cout<<"n\n";
j=-1;
for (i=0; i<numberoffields; i++)
 if (selected[i] && overallbest.weights[++j])
  cout<<overallbest.weights[j]<<endl;
  //cout<<overallbest.weights[j]/(double)maxintweight<<endl;
cout<<"\nTraining file results (migration number/space/total penalty):\n";
for (i=0; i<results; i++)
 if (newbestarray[i])
  cout<<endl<<i+1<<" "<<-1 * trainingpenalty[i];
cout<<endl;
cout<<"\nValidation file results (migration number/space/total penalty):\n";
for (i=0; i<results; i++)
 if (newbestarray[i])
  cout<<endl<<i+1<<" "<<-1 * validationpenalty[i];
cout<<endl;
for (i=0; i<results; i++)
 cout<<"\nMigration "<<i+1<<" new overallbest? "
     <<(newbestarray[i] ? "Yes." : "No.");
cout<<endl;
cout<<"\nTraining/Validation pairs:\n";
for (i=0; i<results; i++)
 if (newbestarray[i])
  cout<<-1*trainingpenalty[i]<<" "<<-1*validationpenalty[i]<<endl;
//endtime = clock();
endtime = gettimeofday_dbl();
cout<<"\nElapsed time = "<<endtime - starttime<<endl;

remove_popfiles();
remove_hashfiles();
destroy_sema();

return 0;
}


string prediction(double test_cases[windowsize][MAXFEATURES],
                  histostruct * histolist,
                  int numberoffields,
                  int weights[MAXFEATURES],
                  double mins[MAXFEATURES],
                  double maxes[MAXFEATURES],
                  double ranges[MAXFEATURES])
{
 struct typeinfo
  {
   string name;
   double probability;
   typeinfo * next;
  };
 typeinfo * labels[windowsize];
 typeinfo * typeptr;
 typeinfo * before;
 typeinfo * after;
 double threshold = 0.02; //for no match
 double max_likelihood;
 double current_likelihood;
 double proportion;
 int i, j;
 int index;
 int maxvoteindex;
 string best_type;
 histostruct * aux;
//cout<<"\nWeights = ";
//for (j=0; j<numberoffields; j++)
 //cout<<weights[j]<<" ";
//cout<<endl;
//cout<<"\nTest case "<<testcasenumber<<":\n";
//for (i=0; i<windowsize; i++)
 //{
  //for (j=0; j<numberoffields; j++)
   //cout<<test_cases[i][j]<<" ";
  //cout<<endl;
 //}
for (i=0; i<windowsize; i++)
 {
 labels[i] = 0;
 aux = histolist;
 while (aux != 0)
  {
   current_likelihood = 1.0;
   //current_likelihood = 0;
   for (j=0; j<numberoffields; j++)
    {
     if (!weights[j]) continue; //most will be 0
     if (test_cases[i][j] > maxes[j])
       proportion = 1.0;
     else if (test_cases[i][j] < mins[j])
       proportion = 0.0;
     else if (maxes[j] == mins[j])
       proportion = 1.0; 
     else
       proportion = (test_cases[i][j]-mins[j])/
                        ranges[j];
     if (proportion >= (1.0*BINS-1.0)/BINS)
        h_index = BINS - 1;
       else
        h_index = ((int)(proportion * BINS)) % BINS;
      current_likelihood *= aux -> histogram[j][h_index];
      //current_likelihood *= aux -> histogram[j][h_index] * weights[j];
     //current_likelihood *= pow(aux -> histogram[j][h_index], (double)weights[j]/maxintweight);
     //current_likelihood += log(aux -> histogram[j][h_index]) * (double)weights[j]/maxintweight;
    }
//insert new type in descending order by probability
 typeptr = new typeinfo;
 typeptr -> name = aux -> type;
 typeptr -> probability = current_likelihood;
 if (labels[i] == 0)
  {
   typeptr -> next = 0;
   labels[i] = typeptr;
  }
 else
  {
   before = labels[i];
   after = 0;
   while (before && (before->probability > current_likelihood))
    {
     after = before;
     before = before->next;
    }
   if (!after)
    {
     typeptr->next = before;
     labels[i] = typeptr;
    }
   else
    {
     typeptr->next = before;
     after->next = typeptr;
    }
  }
 aux = aux -> next;
   }
}
struct votestruct
{
 string type;
 double points;
};
votestruct votes[MAXTYPES];
bool typefound;
int typesfound = 0;
int points, pos;
string currentlabel;
 for (i=0; i<windowsize; i++)
  {
   typeptr = labels[i];
   for (pos=0; pos<maxpoints; pos++)
    {
     currentlabel = typeptr->name;
     typefound = 0;
     for (j=0; j<typesfound; j++)
      if (votes[j].type == currentlabel)
       {
        typefound = 1;
        votes[j].points += maxpoints - pos;
        break;
       }
       if (!typefound)
        {
         votes[typesfound].type = currentlabel;
         votes[typesfound].points = maxpoints - pos;
         typesfound++;
        }
      //for (int jj=0; jj<typesfound; jj++)
       //cout<<"\nvotes["<<jj<<"].type = "<<votes[jj].type
           //<<", points = "<<votes[jj].points;
      typeptr = typeptr -> next;
     }
   }

 maxvoteindex = 0;
 for (i=1; i<typesfound; i++)
  if (votes[i].points > votes[maxvoteindex].points)
   maxvoteindex = i;
 best_type = votes[maxvoteindex].type;
for (i=0; i<windowsize; i++)
 {
  after = labels[i];
  while (labels[i])
   {
    labels[i] = labels[i]->next;
    delete after;
    after = labels[i];
   }
 }
//cout<<"\nbest_type = "<<best_type;
 return best_type;
}

void updatehistos(double examplebuffer[],
             string  class_label,
             int numberoffields,
             double ranges[],
             double min_values[],
             double max_values[],
             histostruct * &histolist)
{
//cout<<"\nIncoming class_label: "<<class_label;
//cout<<"\nfield values: \n";
//for (int z=0; z<numberoffields; z++)
 //cout<<examplebuffer[z]<<", ";
//cout<<flush;
 int i, j, k;
 double proportion;
 int index;
 histostruct * aux;
 bool found = false;
 aux = histolist;
 while (aux != 0)
  {
   if (class_label == aux -> type)
    {
     found = true;
     break;
    }
   aux = aux -> next;
  }
 if (!found)
   {
    aux = new histostruct;
    aux -> type = class_label;
    aux -> next = histolist;
    histolist = aux;
    for (j=0; j<numberoffields; j++)
     for (k=0; k<BINS; k++)
      //aux -> histogram[j][k] = 0.125;
      //aux -> histogram[j][k] = 0.03125;
      aux -> histogram[j][k] = 0;
   }
  //else
    for (j=0; j<numberoffields; j++)
      {
       if (ranges[j] == 0.0)
        proportion = 1.0;
       else
        proportion = (examplebuffer[j]-min_values[j])/ranges[j];
       if (proportion >= (1.0*BINS-1.0)/BINS)
        index = BINS - 1;
       else if (proportion < 0)
        index = 0;
       else
        index = ((int)(proportion * BINS)) % BINS;
       aux -> histogram[j][index]++;
if (index<0 || index > BINS-1)
{
cout<<"\nindex = "<<index;
cout<<"\nproportion = "<<proportion;
cout<<"\nmin_values["<<j<<"] = "<<min_values[j];
cout<<"\nmax_values["<<j<<"] = "<<max_values[j];
cout<<"\nranges["<<j<<"] = "<<ranges[j];
cout<<"\nexamplebuffer["<<j<<"] = "<<examplebuffer[j];
}
      }
}

void rank(int popsize, int numbest)
{
int g,h;
int top_rank;
double max_fitness;
int minindex;
chromosome tempstruct; //for sorting
//sort population array
for (g=0;g<popsize-1;g++)
 {
  minindex = g;
 for (h=g+1;h<popsize;h++)
  if (population[h].fitness < population[minindex].fitness)
    minindex = h;
  tempstruct = population[g];
  population[g] = population[minindex];
  population[minindex] = tempstruct;
 }
 //update overall best chromosome
 //if (population[popsize-1].fitness > overallbest.fitness)
//{
  //overallbest = population[popsize-1];
//cout<<"\nNew overallbest, fitness = "<<overallbest.fitness<<endl;
//}
 //else if (population[popsize-1].fitness < overallbest.fitness)
  //{
   //population[popsize-1] = overallbest;
//cout<<"\nReplaced current best with overallbest, fitness = "
    //<<overallbest.fitness<<endl;
  //}
//put the numbest best in bestN array
int bestloc = 0;
int poploc = popsize - 1;
int lastloc = popsize - 1;
int toofew = 0;
bestN[bestloc] = population[poploc];
for (g=1; g<numbest; g++)
 {while (equal(poploc,lastloc))
  {
   poploc--;
   if (poploc < 0)
    break;
   }
  if (poploc >= 0)
   {
    //population[poploc] is next best
    lastloc = poploc;
    bestloc++;
    bestN[bestloc] = population[poploc];
   }
  else
   {toofew = 1; break;}
  }
if (toofew)
 for (g=bestloc+1; g<numbest; g++)
  bestN[g] = population[popsize - 1];
  
//rank the population, giving lowest rank to the strings
//with lowest fitness, and giving tied strings same rank
population[0].ranking = 1;
for (g=1;g<popsize;g++)
 if (population[g].fitness == population[g-1].fitness)
 //if (absdiff(population[g].fitness,population[g-1].fitness) < 0.0000001)
                                //is either == or >
  population[g].ranking = population[g-1].ranking;
 else
  population[g].ranking = population[g-1].ranking + 1;

 //calculate total of rank values
 rank_total = 0;
 for (g=0; g< popsize; g++)
  rank_total = rank_total + population[g].ranking;
}

void mate(int gen)
{
//select 2 parents from population[]
//perform 1-point crossover, mutation
//place 2 offspring in newpop[]
 //char aux; //for inversion
 int totalones;
 int temp;
 int i,j,randomnumber;
 int current_offset;
 int point, point1, point2; //for inversion
 int mother = select_chrom();
 int father = select_chrom(); //selection with replacement
 int leftpoint, rightpoint; //for crossover
 int leftdelta, rightdelta, tempint;//for uniform crossover
 int length = numberselected;
 //char tempchar;
 //char newsymbol, currentsymbol;
 int newsymbol, currentsymbol;
 //perform crossover
   if ((rand()%1000)/1000.0 > crossover_rate)
    {
     newpop[new_index+1] = population[mother];
     newpop[new_index+2] = population[father];
    }
   else
    {
     //perform uniform, 1-point, or 2-point crossover
//1-point only
     //randomnumber = (rand()/1000)%3;
     //if (randomnumber == 0)
      //{
       //uniform:
       for (i=0; i<length; i++)
        if ((rand()/1000)%2)
         {
          newpop[new_index+1].weights[i] = population[mother].weights[i];
          newpop[new_index+2].weights[i] = population[father].weights[i];
         }
        else
         {
          newpop[new_index+1].weights[i] = population[father].weights[i];
          newpop[new_index+2].weights[i] = population[mother].weights[i];
         }
      //}
     //else if (randomnumber == 1)
      //{
       //1-point:
       //leftpoint = (rand()/1000) % length;
       //for (i=0; i<=leftpoint; i++)
         //newpop[new_index+1].weights[i] = population[mother].weights[i];
       //for (i=leftpoint+1; i<length; i++)
         //newpop[new_index+1].weights[i] = population[father].weights[i];
       //for (i=0; i<=leftpoint; i++)
         //newpop[new_index+2].weights[i] = population[father].weights[i];
       //for (i=leftpoint+1; i<length; i++)
         //newpop[new_index+2].weights[i] = population[mother].weights[i];
       for (j=new_index+1; j<=new_index+2; j++)
        {
         totalones = 0;
         for (i=0; i<length; i++)
          totalones += newpop[j].weights[i];
         if (totalones > subsetsize) //change ones to zeros
          while (totalones > subsetsize)
           {
            point = (rand()/100)%length;
            while (newpop[j].weights[point] == 0)
             point = (rand()/100)%length;
            newpop[j].weights[point] = 0;
            totalones--;
           }
         else if (totalones < subsetsize) //change zeros to ones
          while (totalones < subsetsize)
           {
            point = (rand()/100)%length;
            while (newpop[j].weights[point] == 1)
             point = (rand()/100)%length;
            newpop[j].weights[point] = 1;
            totalones++;
           }
         }
/*
      }
     else //if (randomnumber == 2)
      {
       //2-point:
       leftpoint = (rand()/1000) % (length-1);//max=length-2
       rightpoint =
        ((rand()/1000) % ((length-leftpoint)-1)) + leftpoint + 1;
      for (i=0; i<=leftpoint; i++)
       {
        newpop[new_index+1].string[i] = population[mother].string[i];
        newpop[new_index+2].string[i] = population[father].string[i];
       }
      for (i=leftpoint+1; i<rightpoint; i++)
       {
        newpop[new_index+1].string[i] = population[father].string[i];
        newpop[new_index+2].string[i] = population[mother].string[i];
       }
      for (i=rightpoint; i<length; i++)
       {
        newpop[new_index+1].string[i] = population[mother].string[i];
        newpop[new_index+2].string[i] = population[father].string[i];
       }
     }
*/
    }
 //mutate strings
//swap 2 to preserve subset size
  for (j=new_index+1; j<=new_index+2; j++)
   {
      if ((rand()%1000)/1000.0 <= mutation_rate)
       {
        point1 = (rand()/100)%length;
        point2 = (rand()/100)%length;
        //while (point1 == point2)
        while (newpop[j].weights[point1] ==
               newpop[j].weights[point2])
         point2 = (rand()/100)%length;
        temp = newpop[j].weights[point1];
        newpop[j].weights[point1] = newpop[j].weights[point2];
        newpop[j].weights[point2] = temp;
       }
  }
 new_index = new_index + 2; //global
}//end mate()

void evaluate(fitnode * fhashtable[buckets], classpenalty * classpenaltylist,
              int member, testcase * testcaselist)
{
 chromosome hypo = population[member];
 int weights[MAXFEATURES];
 int i,j,k;
//cout<<"\nmember = "<<member;
//cout<<"\nweights = ";
//for (i=0; i<numberselected; i++)
 //cout<<hypo.weights[i]<<"_";
 bool agree, found;
for (i=0; i<numberselected; i++)
 weights[i] = hypo.weights[i];
//search fhashtable for locations[]
found = 0;
int bucket = hashfunction(hypo.weights);
fitnode * baux = fhashtable[bucket];
while (baux != 0)
 {
  agree = 1;
  for (j=0; j<numberselected; j++)
    if (baux -> weights[j] != weights[j])
     {
      agree = 0;
      break;
     }
  if (agree == 1)
   {
    found = 1;
    break;
   }
  else
   baux = baux -> next;
 }
//cout<<"\nfound in hash table = "<<found;
if (found == 1)
 population[member].fitness = baux -> fitness;
else
//if not in hash table  -- don't forget to add it to hash table
{
cost = 0.0;
bool found;
testcase * testptr = testcaselist;
//get first windowsize samples
for (int ii=1; ii<windowsize; ii++)
 {
  for (j=0; j<numberselected; j++)
     test_cases[ii][j] = testptr->featurearray[j];
  testptr = testptr->next;
 }
while (testptr)
 {
  //shift down cases
  for (i=0; i<windowsize-1; i++)
   for (j=0; j<numberselected; j++)
    test_cases[i][j] = test_cases[i+1][j];
  for (j=0; j<numberselected; j++)
    test_cases[windowsize-1][j] = testptr->featurearray[j];
  testcasenumber++;
  predicted_label = prediction(test_cases, histolist,
                               numberselected, weights,
                               min_values, max_values, ranges);
      if (predicted_label != testptr->class_label)
        //correct++;
       //correct = correct + reward(classrewardlist, class_label);
{
//lock_sema();
//cout<<"\nBefore call to penalty(), class_label = "<<testptr->class_label;
//cout<<", predicted_label = "<<predicted_label<<endl<<flush;
//unlock_sema();
       cost = cost - penalty(classpenaltylist, testptr->class_label, predicted_label);
}
//else
 //{
  //lock_sema();
  //cout<<"\nclass_label = "<<testptr->class_label;
  //cout<<", predicted_label = "<<predicted_label<<endl;
  //unlock_sema();
 //}
  testptr = testptr->next;
 }
population[member].fitness = cost;
//add to hash table
baux = new fitnode;
for (j=0; j<numberselected; j++)
 baux -> weights[j] = weights[j];
baux -> fitness = cost;
baux -> next = fhashtable[bucket];
fhashtable[bucket] = baux;
// Note that new node.
newhashentry *nhe = new newhashentry;
nhe->bucket = bucket;
nhe->entry = baux;
// Add to list of new nodes.
nhe->next = newchildhashentries;
newchildhashentries = nhe;
}
}

void evaluate_validation(classpenalty * classpenaltylist, chromosome overallbest,
              double validationpenalty[], int results, testcase * validationcaselist)
{
 chromosome hypo = overallbest;
 int weights[MAXFEATURES];
 int i,j,k;
//cout<<"\nmember = "<<member;
//cout<<"\nweights = ";
//for (i=0; i<numberselected; i++)
 //cout<<hypo.weights[i]<<"_";
 bool agree, found;
for (i=0; i<numberselected; i++)
 weights[i] = hypo.weights[i];
double cost = 0;
testcase * testptr = validationcaselist;
//get first windowsize samples
for (int ii=1; ii<windowsize; ii++)
 {
  for (j=0; j<numberselected; j++)
     test_cases[ii][j] = testptr->featurearray[j];
  testptr = testptr->next;
 }
while (testptr)
 {
  //shift down cases
  for (i=0; i<windowsize-1; i++)
   for (j=0; j<numberselected; j++)
    test_cases[i][j] = test_cases[i+1][j];
  j = -1;
  for (j=0; j<numberselected; j++)
    test_cases[windowsize-1][j] = testptr->featurearray[j];
  testcasenumber++;
  predicted_label = prediction(test_cases, histolist,
                               numberselected, weights,
                               min_values, max_values, ranges);
      if (predicted_label != testptr->class_label)
        //correct++;
       //correct = correct + reward(classrewardlist, class_label);
       cost = cost - penalty(classpenaltylist, testptr->class_label, predicted_label);
  testptr = testptr->next;
 }
validationpenalty[results] = cost;
}

double penalty(classpenalty * classpenaltylist,
               string class_label, string predicted_label)
{
 if (!classpenaltylist)
  return 1;
 classpenalty * classlistaux;
 penalties * penaltylistaux;
 classlistaux = classpenaltylist;
 while (classlistaux)
  {
   if (classlistaux->classname == class_label)
    break;
   classlistaux = classlistaux->next;
  }
 if (!classlistaux) //no penalties for this class_label
  return 1;
 penaltylistaux = classlistaux->penaltylist;
 if (penaltylistaux->predictedname == "ALL")
  return penaltylistaux->penalty;
 while (penaltylistaux)
  {
   if (penaltylistaux->predictedname == predicted_label)
    return penaltylistaux->penalty;
   penaltylistaux = penaltylistaux->next;
  }
 return 1;
}

int hashfunction(int weights[])
{
double hashvalue = 0.0;
for (int j=0; j<numberselected; j++)
 hashvalue += weights[j]*j;
return ((int) hashvalue) % buckets;
}

double absdiff(double n, double m)
{
 if (n > m) return n - m;
 else return m - n;
}

int equal(int i, int j)
{
 int yes = 1;
 for (int k=0; k<numberselected; k++)
   if (population[i].weights[k] != population[j].weights[k])
    {
     yes = 0;
     break;
    }
 return yes;
}

int select_chrom()
//return index of chromosome based on vector technique
{
 return select_vector[rand() % rank_total];
}

bool samestring(int a, int b)
{
 int length = numberselected;
 for (int j=0; j<length; j++)
  if (population[a].weights[j]
    !=population[b].weights[j])
   return false;
 return true;
}

void updateglobalhistos(double examplebuffer[],
             int numberoffields,
             double ranges[],
             double min_values[],
             double max_values[],
             double globalhistos[MAXFEATURES][MAXBINS])
{
int i, j, k;
 double proportion;
 int index;
    for (j=0; j<numberoffields; j++)
      {
       if (ranges[j] == 0.0)
        proportion = 1.0;
       else
        proportion = (examplebuffer[j]-min_values[j])/ranges[j];
       if (proportion >= (1.0*BINS-1.0)/BINS)
        index = BINS - 1;
       else
        index = ((int)(proportion * BINS)) % BINS;
       globalhistos[j][index]++;
      }
}

int gauss(double mean, double stddev)
{
return (int) (stddev*gasdev() + mean);
}

float gasdev()
{
static int iset=0;
static float gset;
float fac,rsq,v1,v2;

if (iset == 0){
do {
    v1=2.0*(rand()%1000)/1000.0-1.0;
    v2=2.0*(rand()%1000)/1000.0-1.0;
    rsq=v1*v1+v2*v2;
   } while (rsq >=1.0 || rsq ==0.0);
fac = sqrt(-2.0*log(rsq)/rsq);
gset = v1*fac;
iset = 1;
return v2*fac;
} else{
   iset = 0;
   return gset;
  }
}

int read_popfile(char const * const popfilename,
                 chromosome * const population,
                 int const popsize,
                 int const numselected)
{
  ifstream ipop(popfilename);
  string linestring, fieldstring;
  int weight_value = 0;
  char const *buffer = NULL;
  size_t fieldend = 0;
  int i = -1;

  //popsize = 0;

  getline(ipop, linestring);
  while (!ipop.eof())
   {
    ++i;

    // starts with fitness.
    fieldend = linestring.find_first_of(",",0);
    fieldstring.assign(linestring, 0, fieldend);
    buffer = fieldstring.c_str();
    population[i].fitness = strtod(buffer, NULL);
    linestring.erase(0,fieldend+1);

    // assuming ranking is same as file order; this may be
    // a mis-informed (as it is not at all formally informed)
    // assumption.
    population[i].ranking = i;

    for (int j=0; j<numselected; ++j)
     {
      fieldend = linestring.find_first_of(",",0);
      fieldstring.assign(linestring, 0, fieldend);
      buffer = fieldstring.c_str();
      weight_value = atoi(buffer);
      population[i].weights[j]=weight_value;
      linestring.erase(0,fieldend+1);
     }

    //++popsize;
    getline(ipop, linestring);
   }

  return i;
}

int write_popfile(char const * const popfilename,
                  chromosome const * const population, 
                  int const popsize,
                  int const numselected)
{
  ofstream opop(popfilename, ios::trunc|ios::out);
  int i = 0;

  for (i = 0; i < popsize; ++i)
    {
      opop << population[i].fitness;

      for (int j = 0; j < numselected; ++j)
        opop << ',' << population[i].weights[j];

      opop << endl;
    }

  return i;
}

void generate_popfilenames (void)
{
  int fd = -1;
  char temp_filename[PATH_MAX] = { '\0' };
  unsigned long long now = (unsigned long long)time(NULL);

  for (int c = 0; c < NUMCHILDREN; ++c)
    {
      snprintf(
        temp_filename, sizeof(temp_filename),
        "./pop-tmp_" __FILE__ "-%llx-%04d-XXXXXX",
        now, c);

      fd = mkstemp(temp_filename);
      close(fd);
      fd = -1;

      children_popfiles[c] = temp_filename;
    }
}

void remove_popfiles (void)
{
  for (int c = 0; c < NUMCHILDREN; ++c)
      unlink( children_popfiles[c].c_str() );
}

void append_to_hashtable(char const * const hashtablefilename,
                         int const numselected)
{
  ifstream fin(hashtablefilename);
  string linestring(""), fieldstring;
  int bucket = -1;
  char const *buffer = NULL;
  size_t fieldend = 0;

  if ( !fin.eof() )
    getline(fin, linestring);

  while ( !fin.eof() )
    {
      fitnode *f = new fitnode;

      // First is bucket.
      fieldend = linestring.find_first_of(",",0);
      fieldstring.assign(linestring, 0, fieldend);
      buffer = fieldstring.c_str();
      bucket = atoi(buffer);
      linestring.erase(0, fieldend+1);

      // Then Fitness.
      fieldend = linestring.find_first_of(",",0);
      fieldstring.assign(linestring, 0, fieldend);
      buffer = fieldstring.c_str();
      f->fitness = strtod(buffer,  NULL);
      linestring.erase(0, fieldend+1);

      // And the weights.
      for (int c = 0; c < numselected; ++c)
        {
          fieldend = linestring.find_first_of(",",0);
          fieldstring.assign(linestring, 0, fieldend);
          buffer = fieldstring.c_str();
          f->weights[c] = atoi(buffer);
          linestring.erase(0,fieldend+1);
        }

      // Put the new entry in.
      f->next = fhashtable[bucket];
      fhashtable[bucket] = f;

      // next line
      getline(fin, linestring);
    }

}

void write_new_hashtable_entries(char const * const hashtablefile,
                                 const int numselected)
{
  ofstream fout(hashtablefile, ios::trunc|ios::out);
  newhashentry* &h = newchildhashentries;

  while (h)
  {
    fout << h->bucket << ','
         << h->entry->fitness;

    for (int c = 0; c < numselected; ++c)
      fout << ',' << h->entry->weights[c];

    fout << endl;

    // remove entry from list and delete it
    newhashentry *d = h;
    h = h->next;
    delete d;
  }

}

void generate_hashfilenames (void)
{
  int fd = -1;
  char temp_filename[PATH_MAX] = { '\0' };
  unsigned long long now = (unsigned long long)time(NULL);

  for (int c = 0; c < NUMCHILDREN; ++c)
    {
      snprintf(
        temp_filename, sizeof(temp_filename),
        "./hash-tmp_" __FILE__ "-%llx-%04d-XXXXXX",
        now, c);

      fd = mkstemp(temp_filename);
      close(fd);
      fd = -1;

      children_hashtables[c] = temp_filename;
    }
}

void remove_hashfiles (void)
{
  for (int c = 0; c < NUMCHILDREN; ++c)
      unlink( children_hashtables[c].c_str() );
}

void divvy_population (void)
{
  const int num_chromo_per_child = popsize / NUMCHILDREN;
  const int lastextra = popsize - (num_chromo_per_child * NUMCHILDREN);
  int subpop, c, i = 0;

  cerr << __PRETTY_FUNCTION__ << "  Total/Popsize/Per/Last: "
       << NUMCHILDREN << '/' << popsize << '/'
       << num_chromo_per_child << '/'
       << num_chromo_per_child + lastextra << endl;

  // for each chromo
  for (subpop = 0; subpop < NUMCHILDREN; ++subpop)
    for (c = 0; c < num_chromo_per_child; ++c)
        // Copy all features
        children_populations[subpop][c] = population[i++];
 
  // Copy leftovers
  subpop = NUMCHILDREN - 1;
  for (; i < popsize; ++i, ++c)
    children_populations[subpop][c] = population[i];
}

bool spawn_children(void)
{
   for (int c = 0; c < NUMCHILDREN; ++c)
     {
       pid_t p = fork();


       if (p == 0) // If we are a child.
         {
           setup_child(c);
           return true;
         }

       children_pids[c] = p;
     }

   return false; // we are a parent
}

void setup_child(const int child_num)
{
  int num_chromo_per_child = popsize / NUMCHILDREN;

  if (child_num + 1 == NUMCHILDREN)
    num_chromo_per_child += popsize - (num_chromo_per_child * NUMCHILDREN);

  // Move the population from our child group
  // to our parent group.
  for (int c = 0; c < num_chromo_per_child; ++c)
    population[c] = children_populations[child_num][c];

  popsize = num_chromo_per_child;
  popfilename = children_popfiles[child_num];
  child_number = child_num;
  hashtablefilename = children_hashtables[child_num];

  cerr << "Spawned pid " << getpid() << " popsize = " 
       << popsize << "  file: " << popfilename << endl;
}

void wait_for_children(void)
{
  const int num_chromo_per_child = popsize / NUMCHILDREN;
  const int lastextra = popsize - (num_chromo_per_child * NUMCHILDREN);
  const int max_chromo_per_child = num_chromo_per_child + lastextra;


  for (int c = 0; c < NUMCHILDREN; ++c)
    {
      int status = 0;
      int num_chromo = num_chromo_per_child;

      cerr << "Waiting for pid " << children_pids[c] << endl;
      pid_t res = waitpid(children_pids[c], &status, 0);

      if (c + 1 == NUMCHILDREN)
        num_chromo = max_chromo_per_child;

      cerr << "Reading " << num_chromo << " results from " 
           << children_popfiles[c] << endl;

      read_popfile(
        children_popfiles[c].c_str(),
        &children_populations[c][0],
        num_chromo,
        numberselected);
      append_to_hashtable(
        children_hashtables[c].c_str(),
        numberselected);
    }

}

void combine_children_results (double trainingpenalty[], bool newbestarray[], int results)
{
  const int num_chromo_per_child = popsize / NUMCHILDREN;
  const int lastextra = popsize - (num_chromo_per_child * NUMCHILDREN);
  const int max_chromo_per_child = num_chromo_per_child + lastextra;
  int migrationplan[NUMCHILDREN];
  bool used[NUMCHILDREN];
  bool newoverallbest = false;
  int destination, i, j, popnumber;

  cerr << __PRETTY_FUNCTION__ << endl;
  if (results)
   previousbest = overallbest;

  for (int c = 0; c < NUMCHILDREN; ++c)
    {
      int num_chromo = num_chromo_per_child;

      if (c + 1 == NUMCHILDREN)
        num_chromo = max_chromo_per_child;

      // Re-integrate into global population.
      int i = num_chromo_per_child * c;
      for (int j = 0; j < num_chromo; ++i, ++j)
        population[i] = children_populations[c][j];
    }
  for (j=0; j<NUMCHILDREN; j++)
   used[j] = false;
  for (j=0; j<NUMCHILDREN-1; j++)
   {
    popnumber = (rand()/100)%NUMCHILDREN;
    while (popnumber==j || used[popnumber])
     popnumber = (rand()/100)%NUMCHILDREN;
    used[popnumber] = true;
    migrationplan[j] = popnumber;
   }
  for (j=0; j<NUMCHILDREN; j++)
   if (!used[j])
    {
     migrationplan[NUMCHILDREN-1] = j;
     break;
    }
  //insert best member of children_populations[j] at position numbest
  //of children_populations[migrationplan[j]]
  //use numbest in case we keep/restore "elitism"
  for (j=0; j<NUMCHILDREN; j++)
   {
    destination = num_chromo_per_child * migrationplan[j] + numbest;
    population[destination] = children_populations[j][num_chromo_per_child-1];
    if (population[destination].fitness > overallbest.fitness)
     {
      overallbest = population[destination];
      cout<<"\nNew overallbest, fitness = "<<overallbest.fitness<<endl;
      newoverallbest = true;
     }
    cout<<"transferred string from pop "<<j+1<<" to pop "<<migrationplan[j]+1<<" =\n";
      for (i=0; i<numberselected; i++)
       cout<<population[destination].weights[i]<<"_";
      cout<<", "<<population[destination].fitness<<endl;
   }
  newbestarray[results] = newoverallbest;
  trainingpenalty[results] = overallbest.fitness;

}

double gettimeofday_dbl(void)
{
   struct timeval now;
   double nowd;

   gettimeofday(&now, NULL);

   nowd  = now.tv_usec / 1000000.0;
   nowd += now.tv_sec;

   return nowd;
}

int create_sema(void)
{
  // Create new semaphore set with one semaphore.
  semid = semget(IPC_PRIVATE, 1,
    IPC_CREAT | IPC_EXCL |
    S_IRWXU); // Read + Write (+ no-op Excute)
  if (semid < 0)
    return -1;

  // Initialize to 1;
  union semun {
    // I guess the man page says I have to declare this.
    // Which is a little odd...
    int val;
    struct semid_ds *buf;
    unsigned short *array;
  } arg;
  arg.val = 1;
  int res = semctl(semid, 0, SETVAL, arg);

  if (res == 0)
    return 0;

  int local_errno = errno;
  destroy_sema();
  errno = local_errno;
  return -1;
}

void destroy_sema(void)
{
  if (semid < 0)
    return;

  semctl(semid, 0, IPC_RMID);
  semid = -1;
}

inline int mod_sema(int m)
{
  struct sembuf sop;

  sop.sem_num = 0;
  sop.sem_op  = m;
  sop.sem_flg = SEM_UNDO;

  return semop(semid, &sop, 1);
}

void lock_sema(void)
{
  mod_sema(-1);
}

void unlock_sema(void)
{
  mod_sema(+1);
}
