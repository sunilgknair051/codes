#include "Grid.hpp"

#include <omp.h>

Grid::Grid(const u_int n1, const u_int n2):nx_(n1), ny_(n2), sx_(n1 +1), sy_(n2 + 1)
{

	const u_int Lim = sx_ * sy_;

	(*this).data_ = new num_type[Lim];

	#pragma omp parallel for schedule(static)
	for(u_int j =0; j<sy_; ++j)
	{
		for(u_int i =0; i<sx_; ++i)
		{
			(*this).data_[j*sx_ + i] = 0.0;	
		}
	}
}


Grid::~Grid()
{
	delete[] data_;

}

num_type Grid::operator()(const u_int i, const u_int j) const
{
//	assert((i>=0 && j>=0) ||"Grid indices can't be negative");
//	assert((i<sx_ && j<sy_) ||"Grid indices out of range");
	return this->data_[j*sx_ + i];
}


num_type& Grid::operator()(const u_int i, const u_int j)
{
//	assert((i>=0 && j>=0) ||"Grid indices can't be negative");
//	assert((i<sx_ && j<sy_) ||"Grid indices out of range");
	return this->data_[j*sx_ + i];
}

void Grid::computeRHS()
{
	num_type x,y;
	const num_type hx = 2.0/this->nx_;
	const num_type hy = 1.0/this->ny_;
	const u_int nx = this->nx_;
	const u_int ny = this->ny_; 
	
	for(u_int j =1; j<ny; ++j)
	{
		y = j*hy;
		for(u_int i =1; i<nx; ++i)
		{
			x = i*hx;
			(*this)(i,j) = 4*M_PI*M_PI*sin(2*M_PI*x) * sinh(2*M_PI*y);
		}
	}
}

void Grid::applyDirichletBC()
{

	num_type x;
	const num_type hx = 2.0/this->nx_;
	const u_int sx = this->sx_;
	const u_int sy = this->sy_; 
	
	// now aaplying Dirichlet BC on the top horizontal wall
	for(u_int i=0;i<sx;++i)
	{
		x = i*hx;
		(*this)(i,sy-1) = sin(2*M_PI*x) * sinh(2*M_PI);
	}
}


void Grid::gsPoissonSolver(const u_int n_iter, const Grid& f)
{
	const num_type hx = 2.0/this->nx_;
	const num_type hy = 1.0/this->ny_;
	const u_int iLim = this->nx_;
	const u_int jLim = this->ny_; 
	const num_type b = 1.0/(hx*hx);
	const num_type c = 1.0/(hy*hy);
	const num_type a = (k*k + 2*b + 2*c);
	const num_type a_inv = 1.0/a; 
	for(u_int ctr=0; ctr<n_iter; ++ctr)
	{
		for(u_int j=1;j<jLim;++j)
		{
			for(u_int i=1; i<iLim; ++i)
			{
				(*this)(i,j) = a_inv*(f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)));
			}
		}
	}
}



void Grid::rbgsPoissonSolver(const u_int n_iter, const Grid& f)
{
	const num_type hx = 2.0/this->nx_;
	const num_type hy = 1.0/this->ny_;
	const u_int iLim = this->nx_;
	const u_int jLim = this->ny_; 
	const num_type b = 1.0/(hx*hx);
	const num_type c = 1.0/(hy*hy);
	const num_type a = (k*k + 2*b + 2*c);
	const num_type a_inv = 1.0/a;
//	bool jflag = 0;
//	bool bflag = 0;
//	bool iflag = 0;
//	const u_int block = 682;
//	u_int block_num;
//	u_int i_end;
	bool flag =0;
	


	for(u_int ctr=0; ctr<n_iter; ++ctr)
	{

		////////// using flag without if condition : OPT-2
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	


	///////////using if condition for flag : OPT-1
///////////////////////////////////////////////////////////////////////
/*
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
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	}

}

num_type Grid::computeResError(const Grid& f)
{
	const num_type hx = 2.0/this->nx_;
	const num_type hy = 1.0/this->ny_;
	const u_int iLim = this->nx_;
	const u_int jLim = this->ny_; 
	const u_int Lim = (iLim-1) * (jLim-1);
	num_type err = 0, errNorm = 0;
	const num_type b = 1.0/(hx*hx);
	const num_type c = 1.0/(hy*hy);
	const num_type a = (k*k + 2*b + 2*c);
	
	for(u_int j=1;j<jLim;++j)
		{
			for(u_int i=1; i<iLim; ++i)
			{
				err = f(i,j) + b*((*this)(i-1,j)+(*this)(i+1,j)) + c*((*this)(i,j-1)+(*this)(i,j+1)) - a*(*this)(i,j);
				errNorm += err * err; 
			}
		}
		return std::sqrt(errNorm/Lim);

}

void Grid::computeExactSol()
{
	const u_int sx = this->sx_;
	const u_int sy = this->sy_;
	const num_type hx = 2.0/this->nx_;
	const num_type hy = 1.0/this->ny_;
	num_type x,y;
	
	for(u_int j=0; j<sy; ++j)
	{
		y = j*hy;
		for(u_int i=0; i<sx; ++i)
		{
			x = i*hx;
			(*this)(i,j) = sin(2*M_PI*x) * sinh(2*M_PI*y);
		}
	}

}

void Grid::writeGrid(const std::string& filename)
{
	std::ofstream file_writer;
	file_writer.open(filename);
	const u_int sx = this->sx_;
	const u_int sy = this->sy_; 
	const num_type hx = 2.0/ this->nx_;
	const num_type hy = 1.0/ this->ny_;
	num_type x,y;
	
	file_writer<<"# x y u(x,y)"<<std::endl;
	for(u_int j=0; j<sy;++j)
	{
		y = j*hy;
		for(u_int i=0; i<sx;++i)
		{	
			x = i*hx;	
			file_writer<<x<<" "<<y<<" "<<(*this)(i,j)<<std::endl;
		}
		file_writer<<std::endl;
	}
	
	file_writer.close();

}


