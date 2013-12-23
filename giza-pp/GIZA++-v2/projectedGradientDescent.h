#include <algorithm>
#include <vector>
#include <iostream>
//#include "defs.h"
#include <omp.h>

//declaring the global parameters for em with smoothed l0
/*
GLOBAL_PARAMETER(float, ARMIJO_BETA,"armijo_beta","pgd optimization parameter beta used in armijo line search",PARAM_SMOOTHED_LO,0.1);
GLOBAL_PARAMETER(float, ARMIJO_SIGMA,"armijo_sigma","pgd optimization parameter sigma used in armijo line search",PARAM_SMOOTHED_LO,0.0001);
GLOBAL_PARAMETER(float, L0_BETA,"smoothed_l0_beta","optimiztion parameter beta that controls the sharpness of the smoothed l0 prior",PARAM_SMOOTHED_LO,0.5);
GLOBAL_PARAMETER(float, L0_ALPHA,"smoothed_l0_alpha","optimization parameter beta that controls the sharpness of the smoothed l0 prior",PARAM_SMOOTHED_LO,0.0);
*/


using namespace std;


void inline printVector(vector<float> & v)
{
	for (unsigned int i = 0; i< v.size();i++)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}

//template <class PROB>
bool descending (float i,float j) { return (i>j); }

//typedef float PROB ;
//typedef float COUNT ;



void projectOntoSimplex(
    vector<float> &v,
    vector<float> &projection)
{
	vector<float> mu_vector (v.begin(),v.end());
	sort(mu_vector.begin(),mu_vector.end(),descending);
	//vector<float>  projection_step;		
	vector<float> sum_vector (mu_vector.size(),0.0);
	float sum = 0;
	//float max = -99999999;
	vector<float>::iterator it;
	int counter = 1;
	int max_index  = -1;
	float max_index_sum = -99999999;
	for (it = mu_vector.begin() ;it != mu_vector.end(); it++)
	{
		sum += *it;
		float temp_rho = *it - (1/(float)counter)*(sum-1);
		if (temp_rho > 0)
		{
			max_index = counter;
			max_index_sum = sum;
		}
		counter++;
	}
	
	float theta = 1/(float) max_index *(max_index_sum -1);
	//vector <float> projection ;
  
  int index = 0;
	for (it = v.begin() ;it != v.end(); it++)
	{
		float value = max((*it)-theta,(float)0.0);
		projection[index] = value;
    index++;
//		cout<<value<<endl;
	}
	//cout<<"the size of final ector is "<<projection.size()<<endl;
	return;
}

void projectPointsOntoSimplex(
    vector<vector<float> > &new_points,
    vector<vector<float> > &new_feasible_points) {
  int num_conditionals = new_feasible_points.size();
  #pragma omp parallel for firstprivate(num_conditionals) 
  for( int i=0; i<num_conditionals;i++) {
    if (new_points[i].size() > 0) {
      //cout<<"projecting point "<<i<<endl;
      //cout<<"projecting a point with size "<<new_points[i].size()<<endl;
      projectOntoSimplex(new_points[i],new_feasible_points[i]);
    }
  }
  //cout<<"We finished projecting points"<<endl; 
}

float L2_reg_func_value(
       const vector<vector<float> > &  ef_current_points,
       const vector<vector<float> > &  fe_current_points,
       const vector<float> &rowwiseExpCntsSum,
       const vector<vector<pair<unsigned int,unsigned int> > > &ef_map,
       float reg_lambda) {
 	int num_conditionals = ef_current_points.size();
	float func_value = 0.0;
  #pragma omp parallel for firstprivate (num_conditionals) reduction(+ : func_value ) 
	for (int i =0 ;i<num_conditionals;i++)
	{ 
    float conditional_func_value = 0.;
    if (rowwiseExpCntsSum[i] == 0.) {
      continue;
    }
    for (unsigned int j=0; j<ef_map[i].size(); j++) {
      unsigned int f_position =  ef_map[i][j].first;
      unsigned int e_position_in_f_row = ef_map[i][j].second;
      conditional_func_value += (ef_current_points[i][j] - 
                    fe_current_points[f_position][e_position_in_f_row])*
                    (ef_current_points[i][j] - 
                    fe_current_points[f_position][e_position_in_f_row]);
    }
    func_value += conditional_func_value;
  }
  return (func_value*reg_lambda);
   
}
float parallelExpectedLLCompute(
    const vector<vector<float> > & expected_counts,
    const vector<vector<float> > & current_points,
    const vector<float> &rowwiseExpCntsSum) {

 	int num_conditionals = expected_counts.size();
	float func_value = 0.0;

  #pragma omp parallel for firstprivate (num_conditionals) reduction(+ : func_value ) 
	for (int i =0 ;i<num_conditionals;i++)
	{
    // NO NEED TO PERFORM COMPUTATION IF THE EXPECTED COUNTS SUM WAS
    // ZERO
    if (rowwiseExpCntsSum[i] == 0.) {
      continue;
    }
    float conditional_function_value = 0.;
    //cout<<"The size of the row "<<i<<" was "<<current_points[i].size()<<endl;
    for (unsigned int j=0; j<current_points[i].size(); j++) {
      if (current_points[i][j] == 0.0 && expected_counts[i][j] != 0.0)
      {
        cout<<"the probability in position "<<j<<" was 0"<<" and the fractional count was not 0"<<endl;
        exit(0);
      }
      if (current_points[i][j] != 0.0 and expected_counts[i][j] != 0)
      {
        //cout<<"current point is "<<current_points[i][j]<<endl;
        conditional_function_value += expected_counts[i][j] * log(current_points[i][j]);// + L0_ALPHA * exp(-current_point[i]/L0_BETA);
      }
      /*
      else
      {
        func_value += L0_ALPHA * exp(-current_points[i][j]/L0_BETA);
      }
      */
    }
    func_value += conditional_function_value;
	}
	return(func_value);

}

// COMPUTING THE FUNCION VALUE AT THE CURRENT POINT. WE ARE PERFORMING
// PROJECTED GRADIENT DESCENT ON THE NEGATIVE OF THE FUNCTION VALUE
// i.e PROJECTED GRADIENT ASCENT
inline float evalFunction(
    const vector<vector<float> > & ef_expected_counts,
    const vector<vector<float> > & ef_current_points,
    const vector<float> & ef_rowwiseExpCntsSum,
    const vector<vector<float> > & fe_expected_counts,
    const vector<vector<float> > & fe_current_points,
    const vector<float> & fe_rowwiseExpCntsSum,
    const float reg_lambda,
    const vector<vector<pair<unsigned int,unsigned int> > > &ef_map)
{
  //cout<<"In eval function"<<endl;
  float func_value = 0.;

  // FIRST COMPUTING THE EXPECTED COMPLETE DATA LOG LIKELIHOOD IN e given f directiona
  //cout<<"Computing expected complete data LL in the e|f direction"<<endl;
  func_value -= parallelExpectedLLCompute(
    ef_expected_counts,
    ef_current_points,
    ef_rowwiseExpCntsSum);

  //cout<<"Computing expected complete data LL in the f|e direction"<<endl;
  func_value -= parallelExpectedLLCompute(
    fe_expected_counts,
    fe_current_points,
    fe_rowwiseExpCntsSum);
  // Now computing the regularization constant. 
  // For now, the L2 norm
   func_value += L2_reg_func_value(
       ef_current_points,
       fe_current_points,
       ef_rowwiseExpCntsSum,
       ef_map,
       reg_lambda);
	return(func_value);
}

void L2_reg_gradient(
      const vector<vector<float> > &ef_current_points,
      vector<vector<float> > &ef_gradients,
      const vector<vector<float> > &fe_current_points,
      vector<vector<float> > &fe_gradients,
      float reg_lambda,
      const vector<vector<pair<unsigned int,unsigned int> > > &ef_map) {
  
  //First computing the gradient in the ef direction
  int num_conditionals = ef_current_points.size();
  #pragma omp parallel for firstprivate(num_conditionals)
  for(int i=0; i<num_conditionals; i++) {
    for (unsigned int j=0; j<ef_current_points[i].size(); j++) {
      unsigned int f_position =  ef_map[i][j].first;
      unsigned int e_position_in_f_row = ef_map[i][j].second;
      float gradient_value = 2*(ef_current_points[i][j] - 
                    fe_current_points[f_position][e_position_in_f_row]);
      ef_gradients[i][j] += gradient_value;
      fe_gradients[f_position][e_position_in_f_row] -= gradient_value;
    }
  }
}
void evalGradientFromExpLogLikelihood(
    const vector<vector<float> > &expected_counts,
    const vector<float> &rowwiseExpCntsSum,
    const vector<vector<float> > &current_points,
    vector<vector<float> > &gradients) {
  //First evaluating the gradients from the expected complete data log likeliood
	int num_conditionals = expected_counts.size();
  #pragma omp parallel for firstprivate(num_conditionals)
	for (int i =0; i<num_conditionals; i++)
	{
    if (rowwiseExpCntsSum[i] == 0) {
      continue;
    }
    for (unsigned int j=0; j<expected_counts[i].size(); j++) {
      if (current_points[i][j] == 0 && expected_counts[i][j] != 0.0 )
      {
        cout<<"the probability in position "<<j<<" was 0"<<" and the fractional count was not 0"<<endl;
        exit(0);
      }
      if (current_points[i][j] != 0.0)
      {
        gradients[i][j] = -1 *(expected_counts[i][j] / current_points[i][j]) ; //- L0_ALPHA * exp(-current_point[i]/L0_BETA)/L0_BETA	);
      }
    }
  }
}

//template <class COUNT,class PROB>
void inline evalGradient(
      const vector<vector<float> > &ef_expected_counts,
      const vector<vector<float> > &ef_current_points,
      const vector<float> &ef_rowwiseExpCntsSum,
      vector<vector<float> > &ef_gradients,
      const vector<vector<float> > &fe_expected_counts,
      const vector<vector<float> > &fe_current_points,
      const vector<float> &fe_rowwiseExpCntsSum,
      vector<vector<float> > &fe_gradients,
      float reg_lambda,
      const vector<vector<pair<unsigned int,unsigned int> > > &ef_map) {

  evalGradientFromExpLogLikelihood(
      ef_expected_counts,
      ef_rowwiseExpCntsSum,
      ef_current_points,
      ef_gradients);
  evalGradientFromExpLogLikelihood(
      fe_expected_counts,
      fe_rowwiseExpCntsSum,
      fe_current_points,
      fe_gradients);
  L2_reg_gradient(
      ef_current_points,
      ef_gradients,
      fe_current_points,
      fe_gradients,
      reg_lambda,
      ef_map);
}

void zeroInitVectorOfVector(
    const vector<vector<float> > &source,
    vector<vector<float> > &target) {
  int num_conditionals = source.size();
  #pragma omp parallel for firstprivate(num_conditionals)
  for (int i=0; i<num_conditionals; i++) {
    int target_row_size = source[i].size();
    target[i] = vector<float> (target_row_size,0.0);
  }
}

//Initializing the gradients to zero vectors
void initGradientsToZero(
        vector<vector<float> > &ef_gradients,
        const vector<vector<float> > &ef_current_points,
        vector<vector<float> > &fe_gradients,
        const vector<vector<float> > &fe_current_points) {
  zeroInitVectorOfVector(ef_current_points,ef_gradients);
  zeroInitVectorOfVector(fe_current_points,fe_gradients);
}

void getSingleInterpolatedPoints(
    const vector<vector<float> > &new_feasible_points,
    const vector<vector<float> > &current_points,
    vector<vector<float> > &temp_points,
    float current_alpha)  {
  int num_conditionals = new_feasible_points.size();
  #pragma omp parallel for firstprivate(num_conditionals)
  for (int i=0; i<num_conditionals; i++) {
    for (unsigned int j=0; j<new_feasible_points[i].size(); j++) {
      temp_points[i][j] = (1.0 - current_alpha) * current_points[i][j] + current_alpha * new_feasible_points[i][j];
    }
  }

}

void getInterpolatedPoints(
    const vector<vector<float> > &ef_new_feasible_points,
    const vector<vector<float> > &ef_current_points,
    vector<vector<float> > &ef_temp_points,
    const vector<vector<float> > &fe_new_feasible_points,
    const vector<vector<float> > &fe_current_points,
    vector<vector<float> > &fe_temp_points,
    float current_alpha) {

  getSingleInterpolatedPoints(ef_new_feasible_points,
      ef_current_points,
      ef_temp_points,
      current_alpha);

  getSingleInterpolatedPoints(fe_new_feasible_points,
      fe_current_points,
      fe_temp_points,
      current_alpha);

}
void singleGradientStep(
    vector<vector<float> > &new_points,
    const vector<vector<float> > &current_points,
    const vector<vector<float> > &gradients,
    float eta) {
  int num_conditionals = new_points.size();
  #pragma omp parallel for firstprivate(num_conditionals)
  for (int i=0; i<num_conditionals; i++) {
    for (unsigned int j=0; j<new_points[i].size(); j++) {
      new_points[i][j] = current_points[i][j] - eta*gradients[i][j];
    }
  }
}

void takeGradientStep(
    vector<vector<float> > &ef_new_points,
    const vector<vector<float> > &ef_current_points,
    const vector<vector<float> > &ef_gradients,
    vector<vector<float> > &fe_new_points,
    const vector<vector<float> > &fe_current_points,
    const vector<vector<float> > &fe_gradients,
    float eta) {
  //Move the probabilities in the direction of the gradient
  singleGradientStep(
      ef_new_points,
      ef_current_points,
      ef_gradients,
      eta);
  singleGradientStep(
      fe_new_points,
      fe_current_points,
      fe_gradients,
      eta);
}

void projectedGradientDescentWithArmijoRule(
    const vector<vector<float> > & ef_expected_counts,
    const vector<vector<float> > & ef_current_probs,
    const vector<float> & ef_rowwiseExpCntsSum,
    vector<vector<float> > & ef_optimized_probs,
    const vector<vector<float> > & fe_expected_counts,
    const vector<vector<float> > & fe_current_probs,
    const vector<float> & fe_rowwiseExpCntsSum,
    vector<vector<float> > & fe_optimized_probs,
    float reg_lambda,
    const vector<vector<pair<unsigned int,unsigned int> > > &ef_map) {
  float eta = ETA;
  /*
  cout<<"The ef expected counts are "<<ef_expected_counts<<endl;
  cout<<"The fe expected counts are "<<fe_expected_counts<<endl;
  getchar();
  cout<<"The ef current point is "<<ef_current_probs<<endl;
  cout<<"The fe current point is "<<fe_current_probs<<endl;
  getchar();
  */
	//cout<<"projected gradient descent here"<<endl;
  //COPYING THE CURRENT POINT
  vector<vector<float> > ef_current_points(ef_current_probs);
  vector<vector<float> > fe_current_points(fe_current_probs);
  //cout<<"the size of ef current points is "<<ef_current_points<<endl;
  //cout<<"the size of fe current points is "<<fe_current_points<<endl;
	//cout <<"the number of PGD iterations is "<<NUM_PGD_ITERATIONS<<endl;
  cerr<<"Performing PGD iterations"<<endl;
	for (int time = 1; time <= NUM_PGD_ITERATIONS; time++)
	{
    eta /time;
		cout<<"time is"<<time<<endl;
		float current_function_value = evalFunction(
      ef_expected_counts,
      ef_current_points,
      ef_rowwiseExpCntsSum,
      fe_expected_counts,
      fe_current_points,
      fe_rowwiseExpCntsSum,
      reg_lambda,
      ef_map);
    cout<<"The function value was "<<current_function_value<<endl;
    //cerr<<"initializing current gradient"<<endl;
    // INITIAZLIZING THE CURRENT GRADIENTS 
		vector<vector<float> > ef_gradients,fe_gradients;
    ef_gradients = vector<vector<float> >(ef_current_probs.size());
    fe_gradients = vector<vector<float> >(fe_current_probs.size());
    initGradientsToZero(
        ef_gradients,
        ef_current_points,
        fe_gradients,
        fe_current_points);
    //cerr<<"evaluating current gradient"<<endl;
		evalGradient(
      ef_expected_counts,
      ef_current_points,
      ef_rowwiseExpCntsSum,
      ef_gradients,
      fe_expected_counts,
      fe_current_points,
      fe_rowwiseExpCntsSum,
      fe_gradients,
      reg_lambda,
      ef_map);
    
    vector<vector<float> > ef_new_points,fe_new_points;
    ef_new_points = vector<vector<float> >(ef_current_probs.size());
    fe_new_points = vector<vector<float> >(fe_current_probs.size());
    zeroInitVectorOfVector(ef_current_points,
        ef_new_points);
    zeroInitVectorOfVector(fe_current_points,
        fe_new_points);
    //cout<<"The size of ef_new points is "<<ef_new_points<<endl;
    //cout<<"The size of fe_new points is "<<fe_new_points<<endl;
    //cout<<"ef gradients is "<<ef_gradients<<endl;
    //cout<<"fe gradients is "<<fe_gradients<<endl;
    /*
    initGradientsToZero(
        ef_gradients,
        ef_new_points,
        fe_gradients,
        fe_new_points);
    */
		//moving in the opposite direction of the gradient
    //
    takeGradientStep(ef_new_points,
        ef_current_points,
        ef_gradients,
        fe_new_points,
        fe_current_points,
        fe_gradients,
        eta);
    //cout<<"We just took a gradient step"<<endl;
		vector<vector<float> > ef_new_feasible_points,fe_new_feasible_points;
    ef_new_feasible_points = vector<vector<float> >(ef_new_points.size());
    fe_new_feasible_points = vector<vector<float> >(fe_new_points.size());

    zeroInitVectorOfVector(ef_current_points,
        ef_new_feasible_points);
    zeroInitVectorOfVector(fe_current_points,
        fe_new_feasible_points);
    /*
    initGradientsToZero(
        ef_gradients,
        ef_new_feasible_points,
        fe_gradients,
        fe_new_feasible_points);
    */
    //cout<<"Projecting points onto the simplex"<<endl;
		// PROJECTING THE POINTS ON THE SIMPLEX
		projectPointsOntoSimplex(ef_new_points,ef_new_feasible_points);
		projectPointsOntoSimplex(fe_new_points,fe_new_feasible_points);
		//printVector(new_feasible_point);
		//cout<<"feasible point dimension is "<<new_feasible_point.size()<<endl;
		//int num_zero_entries = 0;
		
		//for (int i = 0;i<new_feasible_point.size();i++)
		//{
		//	if (new_feasible_point[i] == 0)
		//	{
		//		num_zero_entries++ ;
		//	}
		//}
    
		//cout<<"the number of zero entries is "<<num_zero_entries<<endl;
		//getchar();
		//cout<<"armijo beta is "<<ARMIJO_BETA<<endl;
		//cout<<"armijo sigma is "<<ARMIJO_SIGMA<<endl;
    /*
    //COMPUTING THE ARMIJO BOUND. FOR NOW, THERE IS NO NEED TO COMPUTE IT. JUST DO LINE SEARCH. 
		float armijo_bound = 0.0;
		for (int i=0 ;i<num_elements;i++)
		{
			float bound_term = ARMIJO_SIGMA * ARMIJO_BETA * gradient[i] * (new_feasible_point[i] - current_point[i]); 	
			//cout<<"the grad is "<<gradient[i]<<" the new feasible point is "<<new_feasible_point[i]<<" current point is "<<current_point[i]<<endl;
			//cout<<"the bound term is "<<bound_term<<endl;
			armijo_bound -= bound_term;
			//cout<<"temp armijo bound "<<armijo_bound<<endl;
		}
		*/
		//getchar();
    bool armijo_bound = 0;
		bool terminate_line_srch = 0;
		int num_steps = 1;
		float current_alpha = ARMIJO_BETA ;
		float final_alpha = 0.0 ; //if the function value does not improve at all, then the armijo beta should be 1
		float current_armijo_bound = armijo_bound;
		float best_func_value = current_function_value;
		bool no_update = 1;
		//cout<<"current function value is "<<current_function_value<<endl;
		//printf ("current function value is %.15f\n",current_function_value);
		while(terminate_line_srch != 1 && num_steps <= 10)
		{	
			//cout<<"current armijo bound is "<<current_armijo_bound<<endl;
			//cout<<"we are in teh while loop"<<endl;
			//cout<<"num steps is "<<num_steps<<endl;
		//	current_beta = 
      vector<vector<float> > ef_temp_points,fe_temp_points;
      ef_temp_points = vector<vector<float> >(ef_new_points.size());
      fe_temp_points = vector<vector<float> >(fe_new_points.size());

      zeroInitVectorOfVector(ef_current_points,
          ef_temp_points);
      zeroInitVectorOfVector(fe_current_points,
          fe_temp_points);

      /*
      initGradientsToZero(
          ef_gradients,
          ef_temp_points,
          fe_gradients,
          fe_temp_points);
      */

      getInterpolatedPoints(ef_new_feasible_points,
        ef_current_points,
        ef_temp_points,
        fe_new_feasible_points,
        fe_current_points,
        fe_temp_points,
        current_alpha);
      //cout<<"After interpolation, the ef point is "<<ef_temp_points<<endl;
      //cout<<"After interpolation, the fe point is "<<fe_temp_points<<endl;

      /*
			for (int i =0 ;i<num_elements;i++)		
			{
				temp_point[i] = (1.0 - current_alpha) * current_point[i] + current_alpha * new_feasible_point[i];
			}
      */
      float func_value_at_temp_point = evalFunction(
        ef_expected_counts,
        ef_temp_points,
        ef_rowwiseExpCntsSum,
        fe_expected_counts,
        fe_temp_points,
        fe_rowwiseExpCntsSum,
        reg_lambda,
        ef_map);
			//cout<<"function value at temp point is "<<func_value_at_temp_point<<endl;
			//printf ("function value at temp point is %.15f and the iteration number is %d \n",func_value_at_temp_point,num_steps);
			//printf ("current alpha is %.15f\n",current_alpha);
			//getchar();
			if (func_value_at_temp_point < best_func_value)
			{
				best_func_value = func_value_at_temp_point;
				final_alpha = current_alpha;
				no_update = 0;
				//cout<<"we arrived at a better function value"<<endl;
				//getchar();
			}
			/*
			if (current_function_value - func_value_at_temp_point >= current_armijo_bound)
			{
				//cout<<"the terminate line src condition was met "<<endl;
				terminate_line_srch = 1;
			}
      */
			current_alpha *= ARMIJO_BETA;
			current_armijo_bound *= ARMIJO_BETA;
			num_steps += 1;
			//getchar();
		}
		//printf ("final alpha was %f\n",final_alpha);
		//cout<<"the value of not update was "<<no_update<<endl;
		//getchar();
	
		//vector<float> next_point ;
		if (no_update == 0)
		{
      //cout<<"Best func value was "<<best_func_value<<endl;
      getInterpolatedPoints(ef_new_feasible_points,
        ef_current_points,
        ef_current_points,
        fe_new_feasible_points,
        fe_current_points,
        fe_current_points,
        final_alpha);
      /*
			//next_point.resize(num_elements);
			for (int i =0 ;i<num_elements;i++)
			{
				float coordinate_point = (1.0 - final_alpha)*current_point[i] + final_alpha * new_feasible_point[i];
				//next_point.push_back(coordinate_point);
				current_point[i] = coordinate_point;
			}
			//current_point = next_point;
      */
		}
		else
		{
			//cout<<" not update was true"<<endl;
			break;
		}

	}		
	//new_prob = current_point;
  //Storing the optimized probs in the new point
  ef_optimized_probs = ef_current_points;
  fe_optimized_probs = fe_current_points;
}



