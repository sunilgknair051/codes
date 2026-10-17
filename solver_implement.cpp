/////////////// using blocking with if condition for flag
//////////////////////////////////////////////////////////////////////
/*
		#pragma omp parallel for schedule(static) private(flag)			
		// now iterating over the red nodes	(even nodes i.e. i+j = even)
		for(u_int ii=1; ii<iLim; ii=ii+block-1)
		{
			i_end = std::min(ii+block, iLim);
			for(u_int j=1; j<jLim; ++j)
			{
				if(j&1)
				{
					for(u_int i=ii; i<i_end; i=i+2)
					{
	//					std::cout<<ii<<"\t"<<i<<"\t"<<j<<"\t"<<end<<std::endl;
						(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
					}
				}		
				else
				{
					for(u_int i=1+ii; i<i_end; i=i+2)
					{
						(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));	
					}
				}
			
			}

		}


	// now iterating over the black nodes	(odd nodes i.e. i+j = odd)
		#pragma omp parallel for schedule(static) private(flag)
		for(u_int ii=1; ii<iLim; ii=ii+block)
		{
			i_end = std::min(ii+block, iLim);
			for(u_int j=1; j<jLim; ++j)
			{
				if(j&1)
				{
					for(u_int i=1+ii; i<i_end; i=i+2)
					{
						(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
					}
				}		
				else
				{
					for(u_int i=ii; i<i_end; i=i+2)
					{
						(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));	
					}
				}
			
			}

		}

*/





		
		////////// blocking without if condition
////////////////////////////////////////////////////////////////////	
/*
	#pragma omp parallel private(jflag,bflag,iflag,block_num,i_end)	
	{							
		std::cout<<"RED"<<std::endl;
		block_num = 0;	
		// now iterating over the red nodes	(even nodes i.e. i+j = even)
		for(u_int ii=1;ii<iLim;ii+=block)
		{	
			++block_num;
//			#pragma omp critical
//			{	
//			std::cout<<"I am thread "<<omp_get_thread_num()<<std::endl;
//			}
			bflag = (block_num&1); // bflag should be = 1 for odd block_num
			i_end = std::min(ii+block,iLim);
			#pragma omp for //schedule(static)
			for(u_int j=1; j<jLim; ++j)	
			{
				jflag = !(j&1) ;   // for j =1, we flag = 0
				iflag = !(bflag^jflag); // iflag should be 0 if both are not same

				for(u_int i=ii+iflag; i<i_end; i=i+2)
				{
					std::cout<<ii<<"\t"<<i<<"\t"<<j<<"\t"<<i_end<<std::endl;
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}		
			}
						
		}

		std::cout<<"BLACK"<<std::endl;
		block_num = 0;
		// now iterating over the black nodes	(odd nodes i.e. i+j = odd)
	
		for(u_int ii=1;ii<iLim;ii+=block)
		{	
			++block_num;

//			#pragma omp critical
//			{	
//			std::cout<<"I am thread "<<omp_get_thread_num()<<std::endl;
//			}			
			bflag = (block_num&1); // bflag should be = 1 for odd block_num
			i_end = std::min(ii+block,iLim);
			#pragma omp for //schedule(static) 
			for(u_int j=1; j<jLim; ++j)	
			{
				jflag = (j&1) ;   // for j =1, we flag = 1
				iflag = !(bflag^jflag);
				for(u_int i=ii+iflag; i<i_end; i=i+2)
				{
					std::cout<<ii<<"\t"<<i<<"\t"<<j<<"\t"<<i_end<<std::endl;
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
		
			}
			
		}
	}
*/





////////// using flag without if condition
////////////////////////////////////////////////////////////////////						

	#pragma omp parallel for schedule(static) private(flag)			
		// now iterating over the red nodes	(even nodes i.e. i+j = even)
		for(u_int j=1; j<jLim; ++j)	
		{
			flag = !(j&1) ;   // for j =1, we flag = 0
				
				for(u_int i=1+flag; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
		
		}


	// now iterating over the black nodes	(odd nodes i.e. i+j = odd)
		#pragma omp parallel for schedule(static) private(flag)
		for(u_int j=1; j<jLim; ++j)	
		{
			flag = (j&1) ;
				
				for(u_int i=1+flag; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
		}
	


	

///////////using if condition for flag
///////////////////////////////////////////////////////////////////////
		#pragma omp parallel for schedule(static)			
		// now iterating over the red nodes	(even nodes i.e. i+j = even)
		for(u_int j=1; j<jLim; ++j)
		{
			if(j&1)
			{
				for(u_int i=1; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
			}		
			else
			{
			for(u_int i=2; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
			}
			
		}


	// now iterating over the black nodes	(odd nodes i.e. i+j = odd)
		#pragma omp parallel for schedule(static)
		for(u_int j=1; j<jLim; ++j)
		{
			if(j&1)
			{
				for(u_int i=2; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
			}		
			else
			{
			for(u_int i=1; i<iLim; i=i+2)
				{
					(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
				}
			}
			
		}


////////////////////////////////////////////////////////////////////
// RBGS SERIAL SOLVER DUMP. Do not delete !!!!!!!

	for(u_int ctr=0; ctr<n_iter; ++ctr)
	{
		// now iterating over the red nodes	(even nodes i.e. i+j = even)
			// red nodes with odd i index 
		for(u_int j=1; j<jLim; j=j+2)
		{
			for(u_int i=1; i<iLim; i=i+2)
			{
				(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
			}
		}
			// red nodes with even i index 
		for(u_int j=2; j<jLim; j=j+2)
		{
			for(u_int i=2; i<iLim; i=i+2)
			{
				(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
			}
		}


	// now iterating over the black nodes	(odd nodes i.e. i+j = odd)
			// red nodes with even i index 
		for(u_int j=1; j<jLim; j=j+2)
		{
			for(u_int i=2; i<iLim; i=i+2)
			{
				(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
			}
		}
			// red nodes with even i index 
		for(u_int j=2; j<jLim; j=j+2)
		{
			for(u_int i=1; i<iLim; i=i+2)
			{
				(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
			}
		}

	}
	

