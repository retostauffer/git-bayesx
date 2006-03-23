
#include "first.h"

#include "ia_ok.h"
#include <algorithm>
#include <iterator>




namespace MCMC
{

	// DEFAULT CONSTRUCTOR:
	IA::IA(void)
	{
		max_ia_order=1;
	}


	// CONSTRUCTOR_1
	IA::IA(const datamatrix & d)
	{
		unsigned i,j;

		data = d; 
		nobs = d.rows();
		nvar = d.cols();

		max_ia_order = 2;
		max_num = nvar*(nvar-1)/2; 

		std::vector<interact> test;
		ia_var = test;

		// ***********compute var_ia ********************************************
		vector <unsigned> help (2);
		for(i=0; i<nvar-1; i++)
		{
			for(j=i+1; j<nvar; j++)
			{
				help[0] = i;
				help[1] = j;

				ia_var.push_back(interact(help));
			}
		}

		occurred = vector<unsigned> (ia_var.size(),0);
	}




	// CONSTRUCTOR_2
	// for interactions of order>2 (some day in future....) 
	IA::IA(unsigned order, const datamatrix & d)
	{
		max_ia_order = order;

		std::vector<interact> test;
		ia_var = test;

		data = d; 
		nobs = d.rows();
		nvar = d.cols();
	}



	// COPY CONSTRUCTOR
	IA::IA(const IA & a)
	{
		nobs = a.nobs;
		nvar = a.nvar;
		data = a.nvar;
		ia_var = a.ia_var;
		occurred = a.occurred;
		
		max_ia_order = a.max_ia_order;			
		max_num = a.max_num;	
	}


	 // OVERLOADED ASSIGNMENT OPERATOR
	const IA & IA::operator=(const IA & a)
	{
		if (this==&a)
		  return *this;

		nobs = a.nobs;
		nvar = a.nvar;
		data = a.nvar;
		ia_var = a.ia_var;
		occurred = a.occurred;
		


		max_ia_order = a.max_ia_order;			
		max_num = a.max_num;	

		return *this;
	}





	

  // FUNCTION: make_ia
  // TASK: creates a new interaction term and adds it to ia_var
	void IA::make_ia (vector<unsigned> terms)
  {
	  unsigned i;
	  unsigned pos;

	  datamatrix var_new (nobs,1,1);
	  double * pvar_new;
	  double * pvar0;
	  double * pvar1;

  
	  pvar_new = var_new.getV();
	  pvar0 = data.getV() + terms[0];
	  pvar1 = data.getV() + terms[1];

	  for(i=0; i<nobs; i++, pvar_new++)
	  {
		  *pvar_new = *pvar0 * *pvar1;

          assert(*pvar_new==0 || *pvar_new==1);

		  pvar0 = pvar0 + nvar;
		  pvar1 = pvar1 + nvar;
	  }

	  pos = get_pos(terms);

	  add_ia(var_new, pos);
  }





  // FUNCTION: choose_ia_term
  // TASK: chooses a new interaction term of order 2
  // which is NOT already in current_ia
	vector<unsigned> IA::choose_ia ( const Matrix<unsigned> & col, 
									 vector <vector <unsigned> > & current_ia)
  {
        unsigned i;
  		bool found;					// help variable
		bool already_there;			// help variable
		unsigned var;				// help variable
		unsigned num;				// number of main effects in the regression model
		vector<unsigned> vec_ia (2);// vector which represents the interaction

		found=false;

        num = 0;
		for(i=0; i<col.rows(); i++)
		{
			if(col(i,0) ==1)
				num++;
		}

		while(found==true)
		{
			// fix first variable of interaction term
			already_there=false;

			while (already_there==true)
			{
				var = rand() % num;

				if(col(var,0)==1)
				{
					vec_ia[0] = var;
					already_there=true;
				}
			}


			// fix secondvariable of interaction term
			already_there=false;

			while (already_there==true)
			{
				var = rand() % num;

				if(var!=vec_ia[0] && col(var,0)==1)
				{
					vec_ia[1]= var;
					already_there=true;
				}
			}
			
			// "sort"
			if(vec_ia[0]>vec_ia[1])
			{
				var = vec_ia[1];
				vec_ia[1] = vec_ia[0];
				vec_ia[0] = var;
			}


			// check if this ia is already in the model
			std:: vector < vector <unsigned > > :: iterator it_i; 
			it_i = current_ia.begin();

			while (*it_i <= vec_ia)
			{
				if(*it_i< vec_ia)
					++it_i;
				else
					found=true;
			}
		} //while

		return vec_ia; 
  }


  




	// FUNCTION: choose_ia_term
	// TASK: chooses a new interaction term of order ord 
	// regardless if it is already in current_ia or not
	vector<unsigned> IA::choose_ia ( const Matrix<unsigned>  & col)
	{
		unsigned i;
		bool already_there;					// help variable
		unsigned var;				// help variable
		unsigned num;				// number of main effects in the regression model
		vector<unsigned> vec_ia (2);	// vector which represents the interaction	 

		num = 0;


		for(i=0; i<col.rows(); i++)
		{
			if(col(i,0) ==1)
				num++;
		}
		
		assert(num>1);


		// fix first variable of interaction term
		already_there=false;

		while (already_there==false)
		{
			var = rand() % nvar;

			if(col(var,0)==1)
			{
				vec_ia[0] = var;
				already_there=true;
			}
		}


		// fix first variable of interaction term
		already_there=false;

		while (already_there==false)
		{
			var = rand() % nvar;;

			if(var!=vec_ia[0] && col(var,0)==1)
			{
				vec_ia[1] = var;
				already_there=true;
			}
		}
		
		// "sort"
		if(vec_ia[0]>vec_ia[1])
		{
			var = vec_ia[1];
			vec_ia[1] = vec_ia[0];
			vec_ia[0] = var;
		}

		return vec_ia; 
	}








   // FUNCTION: already_there
   // TASK: returns true if the interaction vec_ia is already in the current model 	
   // which are represented by current_ia
	bool IA::already_there ( const vector<unsigned> & vec_ia, 
					vector <vector <unsigned> > & current_ia)
	{
		
		bool already_there; 
		already_there=false;
		unsigned size = current_ia.size();

		if(size>0 )
		{
			
			if(current_ia[size-1]>vec_ia)
			{
				unsigned i;
				std::vector <vector <unsigned> > :: iterator it_i;
				
				i=0;
				it_i= current_ia.begin();

				while(i<size)
				{
					if(*it_i<vec_ia)
					{
						i++;
						++it_i;
					}
					else if (*it_i==vec_ia)
					{
						already_there=true;
						i=size;
					}
					else
						i=size;
				}
			}
			else if (current_ia[size-1]==vec_ia)
				already_there=true;
		}


	/*	if(size>0)
		{
			if(current_ia[size-1]==vec_ia)
				already_there=true;
			else if (current_ia[size-1]>vec_ia)
			{
				i=0;
				while (current_ia[i] <= vec_ia)
				{
					if(current_ia[i] < vec_ia)
						i++;
					else
					{
						already_there=true;
						i++;
					}
				}
			}
		}*/
		return already_there;
	}


	// FUNCTION: already_there (vec_ia)
    // TASK: returns true if the interaction vec_ia is already in ia_var 
	bool IA::already_there (const vector<unsigned> & vec_ia)
	{
		if(occurred[get_pos(vec_ia)] ==1)
			return true;
		else 
			return false;
	}











	// FUNCTION: string vec_to_str (terms)
	// TASK: changes terms (vector of numbers) into an ordered (!)string
	ST::string  IA::vec_to_str(vector<unsigned> terms)
	{
		unsigned i;
		ST::string s;

		std::sort(terms.begin(), terms.end());

		for(i=0; i<terms.size(); i++)
			s = s + ST::inttostring(terms[i]);

		return s;
	}




	// FUNCTION: add_ia
	// TASK: adds datamatix ia.ia_dat to the corresponding entry of ia_var
	void IA::add_ia(interact ia) 
	{
		unsigned pos;

		pos = get_pos(ia.ia_term);	
		ia_var[pos].ia_dat = ia.ia_dat;	
		occurred[pos] =1 ;
	}



	
	// FUNCTION: add_ia
	// TASK: adds datamatix ia to interaction at ia_var[pos]
	void IA::add_ia(datamatrix & data, unsigned pos) 
	{
		ia_var[pos].ia_dat = data;
		occurred[pos] = 1 ;
	}


	// FUNCTION: get_pos
	// TASK: gives position of ia if all possible interactions of order 2 
	// of nvar variables are stored in an ordered vector 
	unsigned  IA::get_pos(vector<unsigned> ia)
	{
		
		unsigned k;
		unsigned i = ia[0];
		unsigned j = ia[1];
		unsigned position; 


		assert(ia.size()==2);
		assert(i<j);

		if(i==0)
			position = j-1;
		else
		{
			position = 0; 
			for(k=1; k<i+1; k++)
				position = position + (nvar-k);

			position = position + (j-i)-1;
		}

		return position; 
	}



	// FUNCTION: get_ia
	// TASK: returns pointer to the first element of the matrix of interaction ia
	// regardless if it has already existed before or not
	double * IA::get_ia(vector<unsigned> ia)
	{
		double * p;
		unsigned pos = get_pos(ia);

		if( occurred[pos] !=1)
		{
			make_ia(ia);		// new ia_term is created and added
		}
		p = (ia_var[pos].ia_dat).getV();
		return p;
	}



    // FUNCTION: get_ia
	// TASK: returns k-th element of the matrix of interaction ia
	// regardless if it has already existed before or not
	unsigned IA::get_ia_element(unsigned k, vector<unsigned> ia)
	{
		unsigned ia_elem;
		unsigned pos = get_pos(ia);

		if( occurred[pos] !=1)
		{
			make_ia(ia);		// new ia_term is created and added
		}
		ia_elem = (ia_var[pos].ia_dat)(k,0);
		return int(ia_elem);
	}



		





} //namespace MCMC